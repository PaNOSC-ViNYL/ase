import os
import pickle

from ase.calculators import SinglePointCalculator
from ase.atoms import Atoms
from ase.parallel import rank
from ase.utils import devnull
from ase.neb import NEB


class PickleTrajectory:
    def __init__(self, filename, mode='r', atoms=None, master=None,
                 write_first_image=True):
        self.set_atoms(atoms)
        self.offsets = []
        if master is None:
            master = (rank == 0)
        self.master = master
        self.open(filename, mode)

        if write_first_image and atoms is not None:
            self.write()
        
    def open(self, filename, mode):
        self.fd = filename
        if mode == 'r':
            if isinstance(filename, str):
                self.fd = open(filename, mode + 'b')
            self.read_header()
        elif mode == 'a':
            if isinstance(filename, str):
                self.fd = open(filename, mode + 'b+')
            self.read_header()
        elif mode == 'w':
            if self.master:
                if isinstance(filename, str):
                    if os.path.isfile(filename):
                        os.rename(filename, filename + '.bak')
                    self.fd = open(filename, 'wb')
            else:
                self.fd = devnull
        else:
            raise ValueError('mode must be "r", "w" or "a".')

    def set_atoms(self, atoms=None):
        if atoms is not None and not hasattr(atoms, 'get_positions'):
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

        if isinstance(atoms, NEB):
            neb = atoms
            for image in neb.images:
                self.write(image)
            return

        if len(self.offsets) == 0:
            self.write_header(atoms)

        if atoms.has('momenta'):
            momenta = atoms.get_momenta()
        else:
            momenta = None

        d = {'positions': atoms.get_positions(),
             'cell': atoms.get_cell(),
             'momenta': momenta}


        if atoms.get_calculator() is not None:
            d['energy'] = atoms.get_potential_energy()
            d['forces'] = atoms.get_forces(apply_constraint=False)
            try:
                d['stress'] = atoms.get_stress()
            except NotImplementedError:
                d['stress'] = None

            try:
                if atoms.calc.get_spin_polarized():
                    d['magmoms'] = atoms.get_magnetic_moments()
            except (NotImplementedError, AttributeError):
                pass

        if 'magmoms' not in d and atoms.has('magmoms'):
            d['magmoms'] = atoms.get_initial_magnetic_moments()
            
        if self.master:
            pickle.dump(d, self.fd, protocol=-1)
        self.fd.flush()
        self.offsets.append(self.fd.tell())

    def write_header(self, atoms):
        self.fd.write('PickleTrajectory')
        if atoms.has('tags'):
            tags = atoms.get_tags()
        else:
            tags = None
        d = {'pbc': atoms.get_pbc(),
             'numbers': atoms.get_atomic_numbers(),
             'tags': tags,
             'constraints': atoms.constraints}
        pickle.dump(d, self.fd, protocol=-1)
        self.header_written = True
        self.offsets.append(self.fd.tell())
        
    def close(self):
        self.fd.close()

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
            try:
                magmoms = d['magmoms']
            except KeyError:
                magmoms = None    
            atoms = Atoms(positions=d['positions'],
                          numbers=self.numbers,
                          cell=d['cell'],
                          momenta=d['momenta'],
                          magmoms=magmoms,
                          tags=self.tags,
                          pbc=self.pbc,
                          constraint=[c.copy() for c in self.constraints])
            if 'energy' in d:
                calc = SinglePointCalculator(
                    d['energy'], d['forces'], d['stress'], magmoms, atoms)
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
        return traj[index]
    else:
        return [traj[i] for i in range(len(traj))[index]]

def write_trajectory(filename, images):
    """Write image(s) to trajectory.

    Write also energy, forces, and stress if they are already
    calculated."""

    traj = PickleTrajectory(filename, mode='w')

    if not isinstance(images, (list, tuple)):
        images = [images]
        
    for atoms in images:
        # Avoid potentially expensive calculations:
        calc = atoms.get_calculator()
        if (calc is not None and
            (not hasattr(calc, 'calculation_required') or
             calc.calculation_required(atoms,
                                       ['energy', 'forces', 'stress']))):
            atoms.set_calculator(None)
        
        traj.write(atoms)
        atoms.set_calculator(calc)

    traj.close()
