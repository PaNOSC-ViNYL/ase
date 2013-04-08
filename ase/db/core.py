import os
from datetime import datetime

from ase.atoms import Atoms
from ase.parallel import world
from ase.utils import Lock, OpenLock
from ase.calculators.singlepoint import SinglePointCalculator


def connect(name, type='use_filename_extension', use_lock_file=False):
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

    def write(self, name, atoms, data={}, replace=True):
        if world.rank > 0:
            return

        if name is None:
            name = self.create_random_key(atoms)

        with self.lock:
            self._write(name, atoms, data, replace)

    def _write(self, name, atoms, data, replace):
        pass

    def collect_data(self, atoms):
        #dct = {}#'date': datetime.now()}#, 'user': ..., ...}
        dct = atoms2dict(atoms)
        if atoms.calc is not None:
            dct['calculator'] = {'name': atoms.calc.name,
                                 'parameters': atoms.calc.todict()}
            if len(atoms.calc.check_state(atoms)) == 0:
                dct['results'] = atoms.calc.results
            else:
                dct['results'] = {}
        return dct

    def get(self, names, attach_calculator=False):
        if isinstance(names, str) or names in [0, -1]:
            return self._get([names], attach_calculator)[0]
        if isinstance(names, (list, tuple)):
            return self._get(names, attach_calculator)
        if names == slice(None, None, None):
            return self._get([names], attach_calculator)
        assert 0

    def __getitem__(self, index):
        result = self.get(index)
        if isinstance(result, list):
            return [item[0] for item in result]
        return result[0]


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


def dict2atoms(dct, attach_calculator=False):
    atoms = Atoms(dct['numbers'],
                  dct['positions'],
                  cell=dct['cell'],
                  pbc=dct['pbc'],
                  magmoms=dct.get('magmoms'),
                  )#constraint=...)
    results = dct.get('results')
    if attach_calculator:
        atoms.calc = 117
    elif results:
        atoms.calc = SinglePointCalculator(atoms, **results)

    return atoms
