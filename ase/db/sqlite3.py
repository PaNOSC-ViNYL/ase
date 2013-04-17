from __future__ import absolute_import
import os
import json
import sqlite3
import collections

import numpy as np

from ase.parallel import world
from ase.db import KeyCollisionError
from ase.db.core import NoDatabase, dict2atoms
from ase.db.json import encode, numpyfy


class SQLite3Database(NoDatabase):
    def __init__(self, filename, use_lock_file=False):
        NoDatabase.__init__(self, use_lock_file)
        self.filename = filename
        
    def _write(self, id, atoms, extra, replace):
        conn = sqlite3.connect(self.filename)
        c = conn.execute(
            'select count(*) from sqlite_master where name="systems"')
        if not c.fetchone()[0]:
            conn.execute(
                'create table systems (id text primary key, ' +
                'natoms integer, ' +
                'numbers blob, positions blob, cell blob, pbc integer, ' +
                'initial_magmoms blob, initial_charges blob, ' +
                'constraints text, ' +
                'calculator_name text, calculator_parameters blob, ' +
                'energy real, free_energy real, forces blob, stress blob, ' +
                'magmoms blob, magmom blob, charges blob, charge real, ' +
                'extra text)')
            conn.execute(
                'create table species (atomic_number integer, ' +
                'natoms integer, id text, ' +
                'foreign key (id) references systems (id))')
            conn.execute(
                'create index species_index on species (atomic_number)')
            conn.execute(
                'create table text_keyvals (key text, value text, id text, ' +
                'foreign key (id) references systems (id))')
            conn.execute(
                'create table number_keyvals (key text, value real, id text, ' +
                'foreign key (id) references systems (id))')
            conn.commit()

        if atoms is None:
            row = (id, -1, None, None, None, None, None, None, None, None,
                   None, None, None, None, None, None, None, None, None, None)
        else:
            dct = self.collect_data(atoms)
            row = (id,
                   len(atoms),
                   blob(atoms.numbers),
                   blob(atoms.positions),
                   blob(atoms.cell),
                   blob(atoms.pbc),
                   blob(dct.get('magmoms')),
                   blob(dct.get('charges')),
                   '')
            if 'calculator' in dct:
                row += (dct['calculator']['name'],
                        encode(dct['calculator']['parameters']))
            else:
                row += (None, None)
            if 'results' in dct:
                r = dct['results']
                row += (r.get('energy'),
                        r.get('free_energy'),
                        blob(r.get('forces')),
                        blob(r.get('stress')),
                        blob(r.get('magmoms')),
                        blob(r.get('magmom')),
                        blob(r.get('charges')),
                        r.get('charge'))
            else:
                row += (None, None, None, None, None, None, None, None)
            row += (encode(extra),)
        q = ', '.join('?' * len(row))
        if replace:
            conn.execute('insert or replace into systems values (%s)' % q, row)
        else:
            try:
                conn.execute('insert into systems values (%s)' % q, row)
            except sqlite3.IntegrityError:
                raise KeyCollisionError

        if atoms is not None:
            species = [(int(Z), n, id)
                       for Z, n in collections.Counter(atoms.numbers).items()]
            species.sort()
            conn.executemany('insert into species values (?, ?, ?)', species)

        text_keyvals = []
        number_keyvals = []
        for key, value in extra.items():
            if isinstance(value, str):
                text_keyvals.append([key, value, id])
            elif isinstance(value, (float, int)):
                number_keyvals.append([key, float(value), id])

        if text_keyvals:
            conn.executemany('insert into text_keyvals values (?, ?, ?)',
                             text_keyvals)
        if number_keyvals:
            conn.executemany('insert into number_keyvals values (?, ?, ?)',
                             number_keyvals)

        conn.commit()
        conn.close()
       
    def get_dict(self, id):
        conn = sqlite3.connect(self.filename)
        c = conn.cursor()
        if id in [-1, 0]:
            c.execute('select count(*) from systems')
            assert c.fetchone()[0] == 1
            c.execute('select * from systems')
            row = c.fetchone()
            return [self.row_to_atoms_and_extra(row, attach_calculator)]
        elif names == [slice(None, None, None)]:
            c.execute('select * from systems')
            return [self.row_to_atoms_and_extra(row, attach_calculator)
                    for row in c.fetchall()]
        else:
            systems = []
            for name in names:
                c.execute('select * from systems where name=?', (name,))
                row = c.fetchone()
                systems.append(
                    self.row_to_atoms_and_extra(row, attach_calculator))
            return systems

    def row_to_dict(self, row):
        natoms = row[1]
        dct = {'numbers': deblob(row[2], int),
               'positions': deblob(row[3], shape=(natoms, 3)),
               'cell': deblob(row[4], shape=(3, 3)),
               'pbc': deblob(row[5], bool),
               'magmoms': deblob(row[6]),
               'charges': deblob(row[7]),
               'constraints': ''}
        if row[9] is not None:
            dct['calculator'] = {'name': row[9],
                                 'parameters': numpyfy(json.loads(row[10]))}
            results = {}
            if row[11] is not None:
                results['energy'] = row[11]
            if row[12] is not None:
                results['free_energy'] = row[12]
            if row[13] is not None:
                results['forces'] = deblob(row[13], shape=(natoms, 3))
            if row[14] is not None:
                results['stress'] = deblob(row[14], shape=(3, 3))
            if row[15] is not None:
                results['magmoms'] = deblob(row[15])
            if row[16] is not None:
                results['magmom'] = deblob(row[16])[0]
            if row[17] is not None:
                results['charges'] = deblob(row[17])
            if row[18] is not None:
                results['charge'] = row[18]
            if results:
                dct['results'] = results
        dct['extra'] = numpyfy(json.loads(row[19]))
        return dct


def blob(array):
    """Convert array to blob/buffer object."""

    if array is None:
        return None
    return buffer(array)


def deblob(buf, dtype=float, shape=None):
    """Convert blob/buffer object to ndarray of correct dtype and shape.

    (without creating an extra view)."""

    if buf is None:
        return None
    array = np.frombuffer(buf, dtype)
    if shape is not None:
        array.shape = shape
    return array
