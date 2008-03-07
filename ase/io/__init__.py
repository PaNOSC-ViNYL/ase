from tarfile import is_tarfile
from zipfile import is_zipfile

from ase.atoms import Atoms


def read(filename, index=-1, format=None):
    """Read Atoms object(s) from file.

    Parameters
    ==========
    filename: str
        Name of the file to read from.
    index: int or slice
        If the file contains several configurations, the last configuration
        will be returned by default.  Use index=n to get configuration
        number n (counting from zero).
    format: str
        Used to specify the file-format.  If not given, the file-format
        will be guessed by the `filetype` function.
    """
    p = filename.rfind('@')
    if p != -1:
        try:
            index = string2index(filename[p + 1:])
        except ValueError:
            pass
        else:
            filename = filename[:p]

    if format is None:
        format = filetype(filename)

    if format.startswith('gpw'):
        from gpaw import Calculator
        atoms = Calculator(filename, txt=None).get_atoms()
        atoms.set_calculator(None)
        return atoms
    
    if format == 'xyz':
        from ase.io.xyz import read_xyz
        return read_xyz(filename, index)

    if format == 'traj':
        from ase.io.trajectory import read_trajectory
        return read_trajectory(filename, index)

    if format == 'cube':
        from ase.io.cube import read_cube
        return read_cube(filename, index)

    if format == 'nc':
        from ase.io.netcdf import read_netcdf
        return read_netcdf(filename, index)

    if format == 'gpaw-text':
        from ase.io.gpawtext import read_gpaw_text
        return read_gpaw_text(filename, index)

    if format == 'dacapo-text':
        from ase.io.dacapo import read_dacapo_text
        return read_dacapo_text(filename)

    if format == 'dacapo':
        from ase.io.dacapo import read_dacapo
        return read_dacapo(filename)
    
    raise RuntimeError('That can *not* happen!')


def write(filename, images, format=None, **kwargs):
    """Write Atoms object(s) to file.

    Parameters
    ==========
    filename: str
        Name of the file to write to.
    images: Atoms object or list of Atoms objects
        A single Atoms object or a list of Atoms objects.
    format: str
        Used to specify the file-format.  If not given, the file-format
        will be taken from suffix of the filename.
    """
    
    if format is None:
        if filename is not None:
            suffix = filename.split('.')[-1]
            format = {}.get(suffix, suffix)
        else:
            format = 'xyz'

    if format == 'xyz':
        from ase.io.xyz import write_xyz
        write_xyz(filename, images)
        return

    format = {'traj': 'trajectory', 'nc': 'netcdf'}.get(format, format)
    name = 'write_' + format
    try:
        write = getattr(__import__('ase.io.%s' % format, {}, {}, [name]), name)
    except ImportError:
        raise TypeError('Unknown format: "%s".' % format)
    
    try:
        write(filename, images, **kwargs)
    except TypeError:
        write(filename, images)

def string2index(string):
    if ':' not in string:
        return int(string)
    i = []
    for s in string.split(':'):
        if s == '':
            i.append(None)
        else:
            i.append(int(s))
    i += (3 - len(i)) * [None]
    return slice(*i)


def filetype(filename):
    """Try to guess the type of the file.

    Known formats:

    =========================  ===========
    format                     short name
    =========================  ===========
    GPAW restart-file          gpw
    Dacapo netCDF output file  dacapo
    Old ASE netCDF trajectory  nc
    Virtual Nano Lab file      vnl
    ASE pickle trajectory      traj
    GPAW text output           gpaw-text
    CUBE file                  cube
    Dacapo text output         dacapo-text
    XYZ-file                   xyz
    =========================  ===========
    """
    if is_tarfile(filename):
        return 'gpw-tar'

    fileobj = open(filename)
    if fileobj.read(3) == 'CDF':
        from ase.io.pupynere import NetCDFFile
        nc = NetCDFFile(filename)
        if 'number_of_dynamic_atoms' in nc.dimensions:
            return 'dacapo'

        history = nc.history
        if history == 'GPAW restart file':
            return 'gpw-nc'
        if history == 'ASE trajectory':
            return 'nc'
        if history == 'Dacapo':
            return 'dacapo'
        raise IOError('Unknown netCDF file!')

    if is_zipfile(filename):
        return 'vnl'

    fileobj.seek(0)
    lines = fileobj.readlines(1000)

    if lines[0].startswith('PickleTrajectory'):
        return 'traj'

    if lines[1].startswith('OUTER LOOP:'):
        return 'cube'
    
    if '  ___ ___ ___ _ _ _  \n' in lines:
        return 'gpaw-text'

    if (' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n'
        in lines[:90]):
        return 'dacapo-text'
    return 'xyz'
