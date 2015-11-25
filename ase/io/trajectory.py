from __future__ import print_function
import warnings

from ase.calculators.singlepoint import SinglePointCalculator, all_properties
from ase.constraints import dict2constraint
from ase.atoms import Atoms
from ase.io.aff import affopen, DummyWriter, InvalidAFFError
from ase.io.jsonio import encode, decode
from ase.io.pickletrajectory import PickleTrajectory
from ase.parallel import world

__all__ = ['Trajectory', 'PickleTrajectory']


def Trajectory(filename, mode='r', atoms=None, properties=None, master=None):
    """A Trajectory can be created in read, write or append mode.

    Parameters:

    filename: str
        The name of the file.  Traditionally ends in .traj.
    mode: str
        The mode.  'r' is read mode, the file should already exist, and
        no atoms argument should be specified.
        'w' is write mode.  The atoms argument specifies the Atoms
        object to be written to the file, if not given it must instead
        be given as an argument to the write() method.
        'a' is append mode.  It acts as write mode, except that
        data is appended to a preexisting file.
    atoms: Atoms object
        The Atoms object to be written in write or append mode.
    properties: list of str
        If specified, these calculator properties are saved in the
        trajectory.  If not specified, all supported quantities are
        saved.  Possible values: energy, forces, stress, dipole,
        charges, magmom and magmoms.
    master: bool
        Controls which process does the actual writing. The
        default is that process number 0 does this.  If this
        argument is given, processes where it is True will write.

    The atoms, properties and master arguments are ignores in read mode.
    """
    if mode == 'r':
        return TrajectoryReader(filename)
    return TrajectoryWriter(filename, mode, atoms, properties, master=master)
    
    
class TrajectoryWriter:
    """Writes Atoms objects to a .traj file."""
    def __init__(self, filename, mode='w', atoms=None, properties=None,
                 extra=[], master=None):
        """A Trajectory writer, in write or append mode.

        Parameters:

        filename: str
            The name of the file.  Traditionally ends in .traj.
        mode: str
            The mode.  'r' is read mode, the file should already exist, and
            no atoms argument should be specified.
            'w' is write mode.  The atoms argument specifies the Atoms
            object to be written to the file, if not given it must instead
            be given as an argument to the write() method.
            'a' is append mode.  It acts as write mode, except that
            data is appended to a preexisting file.
        atoms: Atoms object
            The Atoms object to be written in write or append mode.
        properties: list of str
            If specified, these calculator properties are saved in the
            trajectory.  If not specified, all supported quantities are
            saved.  Possible values: energy, forces, stress, dipole,
            charges, magmom and magmoms.
        master: bool
            Controls which process does the actual writing. The
            default is that process number 0 does this.  If this
            argument is given, processes where it is True will write.
        """
        if master is None:
            master = (world.rank == 0)
        self.master = master
        self.atoms = atoms
        self.properties = properties
        
        self.numbers = None
        self.pbc = None
        self.masses = None

        self._open(filename, mode)

    def _open(self, filename, mode):
        if mode not in 'aw':
            raise ValueError('mode must be "w" or "a".')
        if self.master:
            self.backend = affopen(filename, mode, tag='ASE-Trajectory')
            if len(self.backend) > 0:
                r = affopen(filename)
                self.numbers = r.numbers
                self.pbc = r.pbc
        else:
            self.backend = DummyWriter()
                
    def write(self, atoms=None, **kwargs):
        """Write the atoms to the file.

        If the atoms argument is not given, the atoms object specified
        when creating the trajectory object is used.
        
        Use keyword arguments to add extra properties::
            
            writer.write(atoms, energy=117, dipole=[0, 0, 1.0])
        """
        b = self.backend

        if atoms is None:
            atoms = self.atoms

        if hasattr(atoms, 'interpolate'):
            # seems to be a NEB
            neb = atoms
            assert not neb.parallel or world.size == 1
            for image in neb.images:
                self.write(image)
            return
        while hasattr(atoms, 'atoms_for_saving'):
            # Seems to be a Filter or similar, instructing us to
            # save the original atoms.
            atoms = atoms.atoms_for_saving

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
        
        if atoms.has('tags'):
            b.write(tags=atoms.get_tags())

        if atoms.has('momenta'):
            b.write(momenta=atoms.get_momenta())

        if atoms.has('magmoms'):
            b.write(magmoms=atoms.get_initial_magnetic_moments())
            
        if atoms.has('charges'):
            b.write(charges=atoms.get_initial_charges())

        calc = atoms.get_calculator()
        
        if calc is None and len(kwargs) > 0:
            calc = SinglePointCalculator(atoms)

        if calc is not None:
            if not hasattr(calc, 'get_property'):
                calc = OldCalculatorWrapper(calc)
            c = b.child('calculator')
            c.write(name=calc.name)
            for prop in all_properties:
                if prop in kwargs:
                    x = kwargs[prop]
                else:
                    if self.properties is not None:
                        if prop in self.properties:
                            x = calc.get_property(prop, atoms)
                        else:
                            x = None
                    else:
                        try:
                            x = calc.get_property(prop, atoms,
                                                  allow_calculation=False)
                        except (NotImplementedError, KeyError):
                            # KeyError is needed for Jacapo.
                            x = None
                if x is not None:
                    if prop in ['stress', 'dipole']:
                        x = x.tolist()
                    c.write(prop, x)

        info = {}
        for key, value in atoms.info.items():
            try:
                encode(value)
            except TypeError:
                warnings.warn('Skipping "{0}" info.'.format(key))
            else:
                info[key] = value
        if info:
            b.write(info=info)

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
            if all(hasattr(c, 'todict') for c in atoms.constraints):
                b.write(constraints=encode(atoms.constraints))
        if atoms.has('masses'):
            b.write(masses=atoms.get_masses())

    def close(self):
        """Close the trajectory file."""
        self.backend.close()

    def __len__(self):
        return world.sum(len(self.backend))


class TrajectoryReader:
    """Reads Atoms objects from a .traj file."""
    def __init__(self, filename):
        """A Trajectory in read mode.

        The filename traditionally ends in .traj.
        """
        
        self.numbers = None
        self.pbc = None
        self.masses = None

        self._open(filename)

    def _open(self, filename):
        try:
            self.backend = affopen(filename, 'r')
        except InvalidAFFError:
            raise RuntimeError('This is not a valid ASE trajectory file. '
                               'If this is an old-format (version <3.9) '
                               'PickleTrajectory file you can convert it '
                               'with ase.io.trajectory.convert("%s") '
                               'or:\n\n $ python -m ase.io.trajectory %s'
                               % (filename, filename))
        self._read_header()

    def _read_header(self):
        b = self.backend
        if b.get_tag() != 'ASE-Trajectory':
            raise IOError('This is not a trajectory file!')

        if len(b) > 0:
            self.pbc = b.pbc
            self.numbers = b.numbers
            self.masses = b.get('masses')
            self.constraints = b.get('constraints', '[]')

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
                                  for d in decode(self.constraints)],
                      momenta=b.get('momenta'),
                      magmoms=b.get('magmoms'),
                      charges=b.get('charges'),
                      tags=b.get('tags'))
        if 'calculator' in b:
            results = {}
            c = b.calculator
            for prop in all_properties:
                if prop in c:
                    results[prop] = c.get(prop)
            calc = SinglePointCalculator(atoms, **results)
            calc.name = b.calculator.name
            atoms.set_calculator(calc)
        return atoms

    def __len__(self):
        return len(self.backend)

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]


def read_traj(filename, index):
    trj = TrajectoryReader(filename)
    for i in range(*index.indices(len(trj))):
        yield trj[i]


def write_traj(filename, images):
    """Write image(s) to trajectory."""
    trj = TrajectoryWriter(filename, mode='w')
    if isinstance(images, Atoms):
        images = [images]
    for atoms in images:
        trj.write(atoms)
    trj.close()

    
class OldCalculatorWrapper:
    def __init__(self, calc):
        self.calc = calc
        try:
            self.name = calc.name
        except AttributeError:
            self.name = calc.__class__.__name__.lower()
    
    def get_property(self, prop, atoms, allow_calculation=True):
        if not allow_calculation and self.calc.calculation_required(atoms,
                                                                    [prop]):
            return None
        method = 'get_' + {'energy': 'potential_energy',
                           'magmom': 'magnetic_moment',
                           'magmoms': 'magnetic_moments',
                           'dipole': 'dipole_moment'}.get(prop, prop)
        try:
            result = getattr(self.calc, method)(atoms)
        except AttributeError:
            raise NotImplementedError
        return result
        

def convert(name):
    import os
    t = TrajectoryWriter(name + '.new')
    for atoms in PickleTrajectory(name, _warn=False):
        t.write(atoms)
    os.rename(name, name + '.old')
    os.rename(name + '.new', name)
    
    
def main():
    import optparse
    parser = optparse.OptionParser(usage='python -m ase.io.trajectory '
                                   'a1.traj [a2.traj ...]',
                                   description='Convert old trajectory '
                                   'file(s) to new format. '
                                   'The old file is kept as a1.traj.old.')
    opts, args = parser.parse_args()
    for name in args:
        convert(name)
  
        
if __name__ == '__main__':
    main()
