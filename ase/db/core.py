import os
from datetime import datetime

from ase.atoms import Atoms
from ase.parallel import world
from ase.utils import Lock, OpenLock
from ase.calculators.singlepoint import SinglePointCalculator


def connect(name, type='use_filename_extension', use_lock_file=False):
    if type == 'use_filename_extension':
        if name is None:
            type = None
        else:
            type = os.path.splitext(name)[1][1:]

    if type is None:
        return NoDatabase()

    if type == 'json':
        from ase.db.json import JSONDatabase as DB
    elif type == 'sqlite3':
        from ase.db.sqlite3 import SQLite3Database as DB
    else:
        assert 0
    return DB(name, use_lock_file=use_lock_file)


class NoDatabase:
    def __init__(self, use_lock_file=False):
        if use_lock_file:
            self.lock = Lock(name + '.lock')
        else:
            self.lock = OpenLock()

    def write(self, id, atoms, data={}, replace=True):
        if world.rank > 0:
            return

        if id is None:
            id = self.create_random_id(atoms)

        #with self.lock:
        #    self._write(id, atoms, data, replace)
        self.lock.acquire()
        try:
            self._write(id, atoms, data, replace)
        finally:
            self.lock.release()

    def _write(self, id, atoms, data, replace):
        pass

    def collect_data(self, atoms):
        dct = {'timestamp': datetime.now(),
               'username': os.getenv('USER')}
        if atoms is None:
            return dct
        dct.update(atoms2dict(atoms))
        if atoms.calc is not None:
            dct['calculator_name'] = atoms.calc.name
            dct['calculator_parameters'] = atoms.calc.todict()
            if len(atoms.calc.check_state(atoms)) == 0:
                dct['results'] = atoms.calc.results
            else:
                dct['results'] = {}
        return dct

    def get_atoms(self, id, attach_calculator=False, extra=False):
        dct = self.get_dict(id)
        atoms = dict2atoms(dct, attach_calculator)
        if extra:
            atoms.info = dct['extra']
        return atoms

    def __getitem__(self, index):
        if index == slice(None, None, None):
            return [self[0]]
        return self.get_atoms(index)


def atoms2dict(atoms):
    data = {
        'numbers': atoms.numbers,
        'pbc': atoms.pbc,
        'cell': atoms.cell,
        'positions': atoms.positions}
    if atoms.has('magmoms'):
        data['magmoms'] = atoms.get_initial_magnetic_moments()
    if atoms.has('charges'):
        data['charges'] = atoms.get_initial_charges()
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
