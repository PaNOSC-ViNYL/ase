from tarfile import is_tarfile
from zipfile import is_zipfile

from ase.atoms import Atoms
from ase.units import Bohr


def read(filename, index=-1, format=None):
    """Read Atoms object(s) from file.

    filename: str
        Name of the file to read from.
    index: int or slice
        If the file contains several configurations, the last configuration
        will be returned by default.  Use index=n to get configuration
        number n (counting from zero).
    format: str
        Used to specify the file-format.  If not given, the file-format
        will be guessed by the *filetype* function.

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
    XCrySDen Structure File    xsf  
    Dacapo text output         dacapo-text
    XYZ-file                   xyz
    =========================  ===========

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
        import gpaw
        r = gpaw.io.open(filename, 'r')
        positions = r.get('CartesianPositions') * Bohr
        numbers = r.get('AtomicNumbers')
        cell = r.get('UnitCell') * Bohr
        pbc = r.get('BoundaryConditions')
        tags = r.get('Tags')
        magmoms = r.get('MagneticMoments')

        atoms = Atoms(positions=positions,
                      numbers=numbers,
                      cell=cell,
                      pbc=pbc)
        if tags.any():
            atoms.set_tags(tags)
        if magmoms.any():
            atoms.set_magnetic_moments(magmoms)

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
    
    if format == 'xsf':
        from ase.io.xsf import read_xsf
        return read_xsf(filename, index)

    if format == 'vasp':
        from ase.io.vasp import read_vasp
        return read_vasp(filename)
    
    raise RuntimeError('That can *not* happen!')


def write(filename, images, format=None, **kwargs):
    """Write Atoms object(s) to file.

    filename: str
        Name of the file to write to.
    images: Atoms object or list of Atoms objects
        A single Atoms object or a list of Atoms objects.
    format: str
        Used to specify the file-format.  If not given, the file-format
        will be taken from suffix of the filename.

    The accepted output formats:
  
    =========================  ===========
    format                     short name
    =========================  ===========
    ASE pickle trajectory      traj
    CUBE file                  cube
    XYZ-file                   xyz
    Protein Data Bank          pdb
    XCrySDen Structure File    xsf  
    gOpenMol .plt file         plt  
    Python script              py
    Encapsulated Postscript    eps
    Portable Network Graphics  png
    Persistance of Vision      pov
    =========================  ===========
  
    The use of additional keywords is format specific.
  
    The ``cube`` and ``plt`` formats accept (plt requires it) a ``data``
    keyword, which can be used to write a 3D array to the file along
    with the nuclei coordinates.
  
    The ``eps``, ``png``, and ``pov`` formats are all graphics formats,
    and accept the additional keywords::
  
      rotation='', show_unit_cell=0, radii=None, bbox=None, colors=None,
      scale=20
  
    rotation: str
      The rotation angles, e.g. '45x,70y,90z'.
    show_unit_cell: int
      Can be 0, 1, 2 to either not show, show, or show all of the unit cell.
    radii: array / float
      An array of same length as the list of atoms indicating the sphere radii.
      A single float specifies a uniform scaling of the default covalent radii.
    bbox: array
      XXX
    colors: array
      An array of same length as the list of atoms, indicating the rgb color
      code for each atom.
    scale: int
      Number of pixels per Angstrom.
      
    The ``pov`` accepts the additional keywords:
    
    XXX
  
    For ``pov`` the elements of the color array can also be strings, or 4,
    and 5 vectors.
  
    XXX

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
    """Try to guess the type of the file."""
    fileobj = open(filename)
    s3 = fileobj.read(3)
    if len(s3) == 0:
        raise IOError('Empty file: ' + filename)
    
    if is_tarfile(filename):
        return 'gpw'

    if s3 == 'CDF':
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

    for word in ['ANIMSTEPS', 'CRYSTAL', 'SLAB', 'POLYMER', 'MOLECULE']:
        if lines[0].startswith(word):
            return 'xsf'
        
    return 'xyz'
