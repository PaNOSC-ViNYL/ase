from tarfile import is_tarfile
from zipfile import is_zipfile

from ase.atoms import Atoms


def read(filename, index=-1):
    p = filename.rfind('@')
    if p != -1:
        index = string2index(filename[p + 1:])
        filename = filename[:p]
        
    type = filetype(filename)

    if type.startswith('gpw'):
        from gpaw import Calculator
        atoms = Calculator(filename, txt=None).get_atoms()
        atoms.SetCalculator()
        return Atoms(atoms)
    
    if type == 'xyz':
        from ase.io.xyz import read_xyz
        return read_xyz(filename, index)

    if type == 'traj':
        from ase.io.trajectory import read_trajectory
        return read_trajectory(filename, index)

    if type == 'nc':
        from ase.io.netcdf import read_netcdf
        return read_netcdf(filename, index)

    if type == 'gpaw-text':
        from ase.io.gpaw import read_gpaw_text
        return read_gpaw_text(filename, index)

    if type == 'dacapo-text':
        from ase.io.dacapo import read_dacapo_text
        return read_dacapo_text(filename)

    if type == 'dacapo':
        from Dacapo import Calculator
        atoms = Calculator.ReadAtoms(filename)
        atoms.SetCalculator()
        return Atoms(atoms)
    
    raise RuntimeError('That can *not* happen!')


def write(fileobj, images, format=None, **kwargs):
    if isinstance(fileobj, str):
        filename = fileobj
        fileobj = open(fileobj, 'w')
    else:
        filename = None

    if format is None:
        if filename is not None:
            suffix = filename.split('.')[-1]
            format = {}.get(suffix, suffix)
        else:
            format = 'xyz'

    if format == 'xyz':
        from ase.io.xyz import write_xyz
        write_xyz(fileobj, images)
        return

    format = {'traj': 'trajectory', 'nc': 'netcdf'}.get(format, format)
    name = 'write_' + format
    try:
        write = getattr(__import__('ase.io.%s' % format, {}, {}, [name]), name)
    except ImportError:
        raise TypeError('Unknown format: "%s".' % format)
    
    write(filename, images, **kwargs)


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
    if is_tarfile(filename):
        return 'gpw-tar'

    fileobj = open(filename)
    if fileobj.read(3) == 'CDF':
        try:
            from Scientific.IO.NetCDF import NetCDFFile
        except ImportError:
            raise ImportError("This is a netCDF file, but I can't figure "
                              'out what kind without Scientific.IO.NetCDF'
                              'installed')
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
    
    if lines[1] == '  ___ ___ ___ _ _ _  \n':
        return 'gpaw-text'

    if lines[1].startswith('OUTER LOOP:'):
        return 'cube'
    
    if (' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n'
        in lines[:90]):
        return 'dacapo-text'
    return 'xyz'
