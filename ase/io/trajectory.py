from __future__ import print_function
import os

from ase.calculators.singlepoint import SinglePointCalculator, all_properties
from ase.constraints import dict2constraint
from ase.atoms import Atoms
from ase.io.bdf import bdfopen
from ase.io.jsonio import encode
from ase.io.pickletrajectory import PickleTrajectory
from ase.parallel import rank, barrier
from ase.utils import devnull, basestring


class TrajectoryWriter:
    """Writes Atoms objects to a .trj file."""
    def __init__(self, filename, mode='w', atoms=None, properties=None,
                 extra=[], master=None, backup=True):
        """A PickleTrajectory can be created in read, write or append mode.

        Parameters:

        filename:
            The name of the parameter file.  Should end in .traj.

        mode='r':
            The mode.

            'r' is read mode, the file should already exist, and
            no atoms argument should be specified.

            'w' is write mode.  If the file already exists, it is
            renamed by appending .bak to the file name.  The atoms
            argument specifies the Atoms object to be written to the
            file, if not given it must instead be given as an argument
            to the write() method.

            'a' is append mode.  It acts a write mode, except that
            data is appended to a preexisting file.

        atoms=None:
            The Atoms object to be written in write or append mode.

        master=None:
            Controls which process does the actual writing. The
            default is that process number 0 does this.  If this
            argument is given, processes where it is True will write.

        backup=True:
            Use backup=False to disable renaming of an existing file.
        """
        if master is None:
            master = (rank == 0)
        self.master = master
        self.backup = backup
        self.atoms = atoms
        
        self.numbers = None
        self.pbc = None
        self.masses = None

        self._open(filename, mode)

    def _open(self, filename, mode):
        """Opens the file."""
        self.fd = filename
        if mode == 'a':
            exists = True
            if isinstance(filename, basestring):
                exists = os.path.isfile(filename)
                if exists:
                    exists = os.path.getsize(filename) > 0
                if exists:
                    self.fd = open(filename, 'rb')
                    self.read_header()
                    self.fd.close()
                barrier()
                if self.master:
                    self.fd = open(filename, 'ab+')
                else:
                    self.fd = devnull
        elif mode == 'w':
            if self.master:
                if self.backup and os.path.isfile(filename):
                    os.rename(filename, filename + '.old')
                self.backend = bdfopen(filename, 'w', tag='ASE-Trajectory')
        else:
            raise ValueError('mode must be "w" or "a".')

    def write(self, atoms=None, **kwargs):
        """Write the atoms to the file.

        If the atoms argument is not given, the atoms object specified
        when creating the trajectory object is used.
        """
        b = self.backend

        if atoms is None:
            atoms = self.atoms

        if hasattr(atoms, 'interpolate'):
            # seems to be a NEB
            neb = atoms
            assert not neb.parallel
            neb.get_energies_and_forces(all=True)
            for image in neb.images:
                self.write(image)
            return

        if len(b) == 0:
            self.write_header(atoms)
        else:
            if (atoms.pbc != self.pbc).any():
                raise ValueError('Bad periodic boundary conditions!')
            elif len(atoms) != len(self.numbers):
                raise ValueError('Bad number of atoms!')
            elif (atoms.numbers != self.numbers).any():
                raise ValueError('Bad atomic numbers!')

        b.write(positions=atoms.get_positions(),
                cell=atoms.get_cell().tolist())
        
        if atoms.has('momenta'):
            b.write(momenta=atoms.get_momenta())

        if atoms.has('magmoms'):
            b.write(magmoms=atoms.get_initial_magnetic_moments())
            
        if atoms.has('charges'):
            b.write(charges=atoms.get_initial_charges())

        calc = atoms.get_calculator()
        if calc is not None:
            c = b.child('calculator')
            c.write(name=calc.name)
            for p in all:
                if p in kw:

                    x=kw[p]
                elif pr is not None and p in pr:
                    if p not in res:
                        calc.calc(p)
                    x=res[p]
                elif p in res:
                    x=res[p]
                write(p=x)

        #d['info'] = stringnify_info(atoms.info)

        b.sync()
        
    def write_header(self, atoms):
        # Atomic numbers and periodic boundary conditions are only
        # written once - in the header.  Store them here so that we can
        # check that they are the same for all images:
        self.numbers = atoms.get_atomic_numbers()
        self.pbc = atoms.get_pbc()

        b = self.backend
        b.write(version=1,
                pbc=self.pbc.tolist(),
                numbers=self.numbers)
        if atoms.constraints:
            b.write(constraints=encode(atoms.constraints))
        if atoms.has('masses'):
            b.write(masses=atoms.get_masses())

    def close(self):
        """Close the trajectory file."""
        self.backend.close()

    def __len__(self):
        return len(self.backend)


class TrajectoryReader:
    """Reads/writes Atoms objects from/to a .trj file."""
    def __init__(self, filename, properties=None,
                 extra=[], master=None):
        """A PickleTrajectory can be created in read, write or append mode.

        Parameters:

        filename:
            The name of the parameter file.  Should end in .traj.

        mode='r':
            The mode.

            'r' is read mode, the file should already exist, and
            no atoms argument should be specified.

            'w' is write mode.  If the file already exists, it is
            renamed by appending .bak to the file name.  The atoms
            argument specifies the Atoms object to be written to the
            file, if not given it must instead be given as an argument
            to the write() method.

            'a' is append mode.  It acts a write mode, except that
            data is appended to a preexisting file.

        atoms=None:
            The Atoms object to be written in write or append mode.

        master=None:
            Controls which process does the actual writing. The
            default is that process number 0 does this.  If this
            argument is given, processes where it is True will write.

        backup=True:
            Use backup=False to disable renaming of an existing file.
        """
        if master is None:
            master = (rank == 0)
        self.master = master
        
        self.numbers = None
        self.pbc = None
        self.masses = None

        self._open(filename)

    def _open(self, filename):
        self.backend = bdfopen(filename, 'r')
        self._read_header()

    def _read_header(self):
        b = self.backend
        if b.get_tag() != 'ASE-Trajectory':
            raise IOError('This is not a trajectory file!')

        self.pbc = b.pbc
        self.numbers = b.numbers
        self.masses = b.get('masses')
        self.constraints = b.get('constraints', [])

    def close(self):
        """Close the trajectory file."""
        self.backend.close()

    def __getitem__(self, i=-1):
        b = self.backend[i]
        atoms = Atoms(positions=b.positions,
                      numbers=self.numbers,
                      cell=b.cell,
                      masses=self.masses,
                      pbc=self.pbc,
                      info=b.get('info'),
                      constraint=[dict2constraint(d)
                                  for d in self.constraints],
                      momenta=b.get('momenta'),
                      magmoms=b.get('magmoms'),
                      charges=b.get('charges'),
                      tags=b.get('tags'))
        if 'calculator' in b:
            calc = SinglePointCalculator(atoms, **b.calculator.results)
            calc.name = b.calculator.name
            atoms.set_calculator(calc)
        return atoms

    def __len__(self):
        return len(self.backend)

    def __iter__(self):
        return self

    def __next__(self):
        try:
            return self[len(self.offsets) - 1]
        except IndexError:
            raise StopIteration
    
    next = __next__


def read_trajectory(filename, index=-1):
    trj = TrajectoryReader(filename)
    if isinstance(index, int):
        return trj[index]
    else:
        return [trj[i] for i in range(*index.indices(len(trj)))]


def write_trajectory(filename, images):
    """Write image(s) to trajectory."""
    trj = TrajectoryWriter(filename, mode='w')
    if isinstance(images, Atoms):
        images = [images]
    for atoms in images:
        trj.write(atoms)
    trj.close()

    
t = TrajectoryWriter('a.t', 'w')
t.write(Atoms('H'))
t.close()
print(TrajectoryReader('a.t')[0])
