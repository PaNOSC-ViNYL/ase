import os
import pickle

from ase.calculators import SinglePointCalculator
from ase.atoms import Atoms
from ase.parallel import rank
from ase.utils import devnull


class PickleTrajectory:
    def __init__(self, filename, mode='r', atoms=None):
        self.set_atoms(atoms)
        self.offsets = []
        self.open(filename, mode)
        
    def open(self, filename, mode):
        if mode in 'ar':
            self.fd = open(filename, mode + 'b+')
            self.read_header()
        elif mode == 'w':
            if rank == 0:
                if os.path.isfile(filename):
                    os.rename(filename, filename + '.bak')
                self.fd = open(filename, 'wb')
            else:
                self.fd = devnull
        else:
            raise ValueError('mode must be "r", "w" or "a".')

    def set_atoms(self, atoms=None):
        if atoms is not None and not hasattr(atoms, 'get_atomic_numbers'):
            raise TypeError('"atoms" argument is not an Atoms object.')
        self.atoms = atoms

    def read_header(self):
        try:
            if self.fd.read(len('PickleTrajectory')) != 'PickleTrajectory':
                raise IOError('This is not a trajectory file!')
            d = pickle.load(self.fd)
        except EOFError:
            raise EOFError('Bad trajectory file.')
        self.pbc = d['pbc']
        self.numbers = d['numbers']
        self.tags = d.get('tags')
        self.constraints = d['constraints']
        self.offsets.append(self.fd.tell())

    def write(self, atoms=None):
        if atoms is None:
            atoms = self.atoms

        if len(self.offsets) == 0:
            self.write_header(atoms)

        d = {'positions': atoms.get_positions(),
             'cell': atoms.get_cell(),
             'momenta': atoms.get_momenta()}

        if atoms.get_calculator() is not None:
            d['energy'] = atoms.get_potential_energy()
            d['forces'] = atoms.get_forces()
            try:
                d['stress'] = atoms.get_stress()
            except NotImplementedError:
                d['stress'] = None

        if rank == 0:
            pickle.dump(d, self.fd, protocol=-1)
        self.fd.flush()
        self.offsets.append(self.fd.tell())

    def write_header(self, atoms):
        self.fd.write('PickleTrajectory')
        d = {'pbc': atoms.get_pbc(),
             'numbers': atoms.get_atomic_numbers(),
             'tags': atoms.get_tags(),
             'constraints': atoms.constraints}
        pickle.dump(d, self.fd, protocol=-1)
        self.header_written = True
        self.offsets.append(self.fd.tell())
        
    def close(self):
        self.fd.close()

    def __del__(self):
        self.close()
        
    def __getitem__(self, i=-1):
        N = len(self.offsets)
        if 0 <= i < N:
            self.fd.seek(self.offsets[i])
            try:
                d = pickle.load(self.fd)
            except EOFError:
                raise IndexError
            if i == N - 1:
                self.offsets.append(self.fd.tell())
            atoms = Atoms(positions=d['positions'],
                          numbers=self.numbers,
                          cell=d['cell'],
                          momenta=d['momenta'],
                          tags=self.tags,
                          pbc=self.pbc,
                          constraints=[c.copy() for c in self.constraints])
            if 'energy' in d:
                calc = SinglePointCalculator(
                    d['energy'], d['forces'], d['stress'], atoms)
                atoms.set_calculator(calc)
            return atoms

        if i >= N:
            for j in range(N - 1, i + 1):
                atoms = self[j]
            return atoms

        i = len(self) + i
        if i < 0:
            raise IndexError('Trajectory index out of range.')
        return self[i]

    def __len__(self):
        N = len(self.offsets) - 1
        while True:
            self.fd.seek(self.offsets[N])
            try:
                pickle.load(self.fd)
            except EOFError:
                return N
            self.offsets.append(self.fd.tell())
            N += 1

    def __iter__(self):
        del self.offsets[1:]
        return self

    def next(self):
        try:
            return self[len(self.offsets) - 1]
        except IndexError:
            raise StopIteration


def read_trajectory(filename, index=-1):
    traj = PickleTrajectory(filename, mode='r')

    if isinstance(index, int):
        indices = [index]
    else:
        indices = range(len(traj))[index]

    return [traj[i] for i in indices]

def write_trajectory(filename, images):
    traj = PickleTrajectory(filename, mode='w')

    if not isinstance(images, (list, tuple)):
        images = [images]
        
    for atoms in images:
        traj.write(atoms)

    traj.close()
