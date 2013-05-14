from __future__ import absolute_import
import json
import sqlite3

import numpy as np

from ase.db import IdCollisionError
from ase.db.core import NoDatabase, ops
from ase.db.json import encode, numpyfy


init_statements = """\
create table systems (
 id text primary key,
 unique_id text,
 timestamp real,
 username text,
 numbers blob,
 positions blob,
 cell blob,
 pbc integer,
 initial_magmoms blob,
 initial_charges blob,
 masses blob,
 tags blob,
 moments blob,
 constraints text,
 calculator_name text,
 calculator_parameters text,
 energy real,
 free_energy real,
 forces blob,
 stress blob,
 magmoms blob,
 magmom blob,
 charges blob,
 data text); -- contains keywords and key_value_pairs also
create index unique_id_index on systems(unique_id);
create table species (
 Z integer,
 n integer,
 id text,
 foreign key (id) references systems(id));
create index species_index on species(Z);
create table keywords (
 keyword text,
 id text,
 foreign key (id) references systems(id));
create table text_key_values (
 key text,
 value text,
 id text,
 foreign key (id) references systems(id));
create table number_key_values (
 key text,
 value real,
 id text,
 foreign key (id) references systems (id))
"""


class SQLite3Database(NoDatabase):
    def _write(self, id, atoms, keywords, key_value_pairs, data, replace):
        con = sqlite3.connect(self.filename)
        cur = con.execute(
            'select count(*) from sqlite_master where name="systems"')
        if not cur.fetchone()[0]:
            for statement in init_statements.split(';'):
                con.execute(statement)
            con.commit()
                
        if isinstance(atoms, dict):
            dct = atoms
            unique_id = dct['unique_id']
            cur = con.execute('select id from systems where unique_id=?',
                              (unique_id,))
            rows = cur.fetchall()
            if rows:
                id = rows[0][0]
        else:
            dct = self.collect_data(atoms)

        if id is None:
            cur = con.execute('select count(*) from systems')
            nrows = cur.fetchone()[0]
            while id is None:
                id = self.create_random_id(nrows)
                cur = con.execute('select count(*) from systems where id=?',
                                  id)
                if cur.fetchone()[0] == 1:
                    id = None
            
        if atoms is None:
            row = (id, None, None, None, None, None, None, None, None, None,
                   None, None, None, None, None, None, None, None, None, None,
                   None, None, None, None)

        row = (id,
               dct['unique_id'],
               self.timestamp,
               dct['username'],
               blob(dct.get('numbers')),
               blob(dct.get('positions')),
               blob(dct.get('cell')),
               blob(dct.get('pbc')),
               blob(dct.get('magmoms')),
               blob(dct.get('charges')),
               blob(dct.get('masses')),
               blob(dct.get('tags')),
               blob(dct.get('moments')),
               dct.get('constraints'))
        if 'calculator_name' in dct:
            row += (dct['calculator_name'],
                    encode(dct['calculator_parameters']))
        else:
            row += (None, None)
        if 'results' in dct:
            r = dct['results']
            magmom = r.get('magmom')
            if magmom is not None:
                # magmom can be one or three numbers (non-collinear case)
                magmom = np.array(magmom)
            row += (r.get('energy'),
                    r.get('free_energy'),
                    blob(r.get('forces')),
                    blob(r.get('stress')),
                    blob(r.get('magmoms')),
                    blob(magmom),
                    blob(r.get('charges')))
        else:
            row += (None, None, None, None, None, None, None)
        row += (encode({'data': data,
                        'keywords': keywords,
                        'key_value_pairs': key_value_pairs}),)

        q = ', '.join('?' * len(row))
        if replace:
            con.execute('insert or replace into systems values (%s)' % q, row)
        else:
            try:
                con.execute('insert into systems values (%s)' % q, row)
            except sqlite3.IntegrityError:
                raise IdCollisionError

        if atoms is not None:
            count = np.bincount(atoms.numbers)
            unique_numbers = count.nonzero()[0]
            species = [(int(Z), int(count[Z]), id) for Z in unique_numbers]
            con.executemany('insert into species values (?, ?, ?)', species)

        text_key_values = []
        number_key_values = []
        for key, value in key_value_pairs.items():
            if isinstance(value, str):
                text_key_values.append([key, value, id])
            elif isinstance(value, (float, int)):
                number_key_values.append([key, float(value), id])
            else:
                assert 0, value
 
        if text_key_values:
            con.executemany('insert into text_key_values values (?, ?, ?)',
                             text_key_values)
        if number_key_values:
            con.executemany('insert into number_key_values values (?, ?, ?)',
                             number_key_values)
        if keywords:
            con.executemany('insert into keywords values (?, ?)',
                            [(keyword, id) for keyword in keywords])

        con.commit()
        con.close()
       
    def _get_dict(self, id):
        con = sqlite3.connect(self.filename)
        c = con.cursor()
        if id in [-1, 0]:
            c.execute('select count(*) from systems')
            assert c.fetchone()[0] == 1
            c.execute('select * from systems')
        else:
            c.execute('select * from systems where id=?', (id,))
        row = c.fetchone()
        return self.row_to_dict(row)

    def row_to_dict(self, row):
        dct = {'id': row[0],
               'unique_id': row[1],
               'timestamp': row[2],
               'username': row[3],
               'numbers': deblob(row[4], int),
               'positions': deblob(row[5], shape=(-1, 3)),
               'cell': deblob(row[6], shape=(3, 3)),
               'pbc': deblob(row[7], bool)}
        if row[8] is not None:
            dct['magmoms'] = deblob(row[8])
        if row[9] is not None:
            dct['charges'] = deblob(row[9])
        if row[10] is not None:
            dct['masses'] = deblob(row[10])
        if row[11] is not None:
            dct['tags'] = deblob(row[11], int)
        if row[12] is not None:
            dct['moments'] = deblob(row[12], shape=(-1, 3))
        if row[13] is not None:
            dct['constraints'] = numpyfy(json.loads(row[13]))
        if row[14] is not None:
            dct['calculator_name'] = row[14]
            dct['calculator_parameters'] = numpyfy(json.loads(row[15]))
            results = {}
            if row[16] is not None:
                results['energy'] = row[16]
            if row[17] is not None:
                results['free_energy'] = row[17]
            if row[18] is not None:
                results['forces'] = deblob(row[18], shape=(-1, 3))
            if row[19] is not None:
                results['stress'] = deblob(row[19])
            if row[20] is not None:
                results['magmoms'] = deblob(row[20])
            if row[21] is not None:
                results['magmom'] = deblob(row[21])[0]
            if row[22] is not None:
                results['charges'] = deblob(row[22])
            if results:
                dct['results'] = results
            dct.update(numpyfy(json.loads(row[23])))
        return dct

    def _select(self, keywords, cmps, limit, offset,
                explain, verbosity):
        tables = set(['systems'])
        where = []
        if keywords:
            tables.add('keywords')
            for keyword in keywords:
                where.append('systems.id=keywords.id and keywords.keyword=%r' %
                             keyword)
        bad = {}
        for key, op, value in cmps:
            if isinstance(key, int):
                bad[key] = bad.get(key, True) and ops[op](0, value)
        cmps2 = []
        for key, op, value in cmps:
            if key in ['id', 'energy', 'magmom', 'timestamp', 'username',
                       'calculator_name']:
                where.append('systems.%s%s%r' % (key, op, value))
            elif key == 'natoms':
                cmps2.append((key, ops[op], value))
            elif isinstance(key, int):
                if bad[key]:
                    cmps2.append((key, ops[op], value))
                else:
                    tables.add('species')
                    where.append('systems.id=species.id and ' +
                                 'species.Z=%d and species.n%s%d' %
                                 (key, op, value))
            elif isinstance(value, str):
                tables.add('text_key_values')
                where.append('systems.id=text_key_values.id and ' +
                             'text_key_values.key=%r and ' % key +
                             'text_key_values.value%s%r' %
                             (op, value))
            else:
                tables.add('number_key_values')
                where.append('systems.id=number_key_values.id and ' +
                             'number_key_values.key=%r and ' % key +
                             'number_key_values.value%s%r' %
                             (op, value))
                
        sql = 'select systems.* from\n  ' + ', '.join(tables)
        if where:
            sql += '\n  where\n  ' + ' and\n  '.join(where)
        if explain:
            sql = 'explain query plan ' + sql
        if verbosity == 2:
            print(sql)
        con = sqlite3.connect(self.filename)
        cur = con.execute(sql)
        if explain:
            for row in cur.fetchall():
                yield row
        else:
            for row in cur.fetchall():
                if cmps2:
                    numbers = deblob(row[4], int)
                    for key, op, value in cmps2:
                        if key == 'natoms':
                            if not op(len(numbers), value):
                                break
                        elif not op((numbers == key).sum(), value):
                            break
                    else:
                        yield self.row_to_dict(row)
                else:
                    yield self.row_to_dict(row)
        
    def delete(self, ids):
        con = sqlite3.connect(self.filename)
        cur = con.executemany('delete from systems where id=?',
                              ((id,) for id in ids))
        con.commit()
        con.close()


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
