import os
from datetime import datetime

from ase.atoms import Atoms
from ase.utils import Lock, OpenLock
from ase.calculators.singlepoint import SinglePointCalculator


def database(name, type='use_filename_extension', use_lock_file=False):
    if type == 'use_filename_extension':
        type = os.path.splitext(name)[1][1:]

    if type is None:
        return NoDatabase()

    if type == 'json':
        from ase.db.jsondb import JSONDatabase as DB
    elif type == 'sqlite':
        from ase.db.sqlite import SQLiteDatabase as DB
    return DB(name, use_lock_file=use_lock_file)


class NoDatabase:
    def __init__(self, use_lock_file=False):
        if use_lock_file:
            self.lock = Lock(name + '.lock')
        else:
            self.lock = OpenLock()

    def write(self, name, atoms, data={}, overwrite=False):
        with self.lock:
            self._write(name, atoms, data, overwrite)

    def _write(self, name, atoms, data, overwrite):
        pass

    def create_dictionary(self, atoms, data={}):
        dct = {}#'date': datetime.now(), 'user': ..., ...}
        if atoms is not None:
            dct['atoms'] = atoms2dict(atoms)
            if atoms.calc is not None:
                dct['calculator'] = atoms.calc.todict()
                dct['results'] = atoms.calc.results
        for key in data:
            assert key not in dct
        dct.update(data)
        return dct

    def create_atoms_and_data(self, dct, attach_calculator=False):
        atoms = dict2atoms(dct.pop('atoms'))
        results = dct.pop('results', None)
        if attach_calculator:
            atoms.calc = atoms.calc.create_calculator()
        elif results:
            atoms.calc = SinglePointCalculator(atoms, **results)
        return atoms, dct

    def get(self, names, attach_calculator=False):
        if isinstance(names, str) or names is None:
            return self._get([names], attach_calculator)[0]
        else:
            return self._get(names, attach_calculator)


def atoms2dict(atoms):
    data = {
        'numbers': atoms.numbers,
        'pbc': atoms.pbc,
        'cell': atoms.cell,
        'positions': atoms.positions}
    if atoms.has('magmoms'):
        data['magmoms'] = atoms.get_initial_magnetic_moments()
    if atoms.constraints:
        data['constraints'] = repr(atoms.constraints)
    return data


def dict2atoms(data):
    atoms = Atoms(data['numbers'],
                  data['positions'],
                  cell=data['cell'],
                  pbc=data['pbc'],
                  magmoms=data.get('magmoms'),
                  )#constraint=...)
    return atoms
