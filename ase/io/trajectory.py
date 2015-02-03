from __future__ import print_function
import os
import sys
import errno
import pickle
import warnings
import collections


from ase.calculators.singlepoint import SinglePointCalculator, all_properties
from ase.atoms import Atoms
from ase.parallel import rank, barrier
from ase.utils import devnull, basestring

from ase.io.pickletrajectory import PickleTrajectory
import ase.io.bdf as bdf


class TrajectoryWriter:
    """Writes Atoms objects to a .trj file."""
    def __init__(self, filename, mode='w', atoms=None, properties=None,
                 extra=[], master=None, backup=True, backend=bdf):
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
        self.backend = backend
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
                self.backend = self.backend.open(filename, 'wb')
        else:
            raise ValueError('mode must be "r", "w" or "a".')

    def write(self, atoms=None, **kwargs):
        """Write the atoms to the file.

        If the atoms argument is not given, the atoms object specified
        when creating the trajectory object is used.
        """
        b = self.backend

        self._call_observers(self.pre_observers)
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
            elif self.sanitycheck and len(atoms) != len(self.numbers):
                raise ValueError('Bad number of atoms!')
            elif self.sanitycheck and (atoms.numbers != self.numbers).any():
                raise ValueError('Bad atomic numbers!')

        if atoms.has('momenta'):
            momenta = atoms.get_momenta()
        else:
            momenta = None

        d = {'positions': atoms.get_positions(),
             'cell': atoms.get_cell(),
             'momenta': momenta}

        if atoms.get_calculator() is not None:
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

        if 'magmoms' not in d and atoms.has('magmoms'):
            d['magmoms'] = atoms.get_initial_magnetic_moments()
        if 'charges' not in d and atoms.has('charges'):
            charges = atoms.get_initial_charges()
            if (charges != 0).any():
                d['charges'] = charges

        if self.write_info:
            d['info'] = stringnify_info(atoms.info)

        b.sync()
        
        self._call_observers(self.post_observers)

    def write_header(self, atoms):
        # Atomic numbers and periodic boundary conditions are only
        # written once - in the header.  Store them here so that we can
        # check that they are the same for all images:
        self.numbers = atoms.get_atomic_numbers()
        self.pbc = atoms.get_pbc()

        b = self.backend
        b.write(version=1,
                pbc=self.pbc,
                numbers=self.numbers,
                constraints=encode(atoms.constraints))
        if atoms.has('masses'):
            b.write(masses=atoms.get_masses())

    def close(self):
        """Close the trajectory file."""
        self.backend.close()


    def __len__(self):
        return len(self.backend)

    def pre_write_attach(self, function, interval=1, *args, **kwargs):
        """Attach a function to be called before writing begins.

        function: The function or callable object to be called.

        interval: How often the function is called.  Default: every time (1).

        All other arguments are stored, and passed to the function.
        """
        if not isinstance(function, collections.Callable):
            raise ValueError('Callback object must be callable.')
        self.pre_observers.append((function, interval, args, kwargs))

    def post_write_attach(self, function, interval=1, *args, **kwargs):
        """Attach a function to be called after writing ends.

        function: The function or callable object to be called.

        interval: How often the function is called.  Default: every time (1).

        All other arguments are stored, and passed to the function.
        """
        if not isinstance(function, collections.Callable):
            raise ValueError('Callback object must be callable.')
        self.post_observers.append((function, interval, args, kwargs))

    def _call_observers(self, obs):
        """Call pre/post write observers."""
        for function, interval, args, kwargs in obs:
            if self.write_counter % interval == 0:
                function(*args, **kwargs)


class TrajectoryReader:
    """Reads/writes Atoms objects from/to a .trj file."""
    def __init__(self, filename, properties=None,
                 extra=[], master=None, backend=bdf):
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
        self.backend = backend
        
        self.numbers = None
        self.pbc = None
        self.masses = None

        self._open(filename, mode)

    def _open(self, filename, mode):
        """Opens the file."""
        self.backend = self.backend.read(filename, 'r')
        self.read_header()

    def _read_header(self):
        be = self.backend
        if be.get_tag() != 'ASE-TRAJ':
            raise IOError('This is not a trajectory file!')

        self.pbc = be.pbc
        self.numbers = be.numbers
        self.masses = be.get('masses')

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
                      constraint=encode(self.constraints),
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
    trj = Trajectory(filename, mode='r')
    if isinstance(index, int):
        return trj[index]
    else:
        return [trj[i] for i in range(*index.indices(len(trj)))]


def write_trajectory(filename, images):
    """Write image(s) to trajectory."""
    trj = Trajectory(filename, mode='w')
    if isinstance(images, Atoms):
        images = [images]
    for atoms in images:
        trj.write(atoms)
    trj.close()
