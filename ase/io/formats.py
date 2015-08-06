import collections
import functools
import inspect
import os
import sys

"""
format 'abc' abc.py: read_abc, write_abc.  Add to
does_not_accept_a_file_descriptor and stores_multiple_images lists




"""
from ase.atoms import Atoms
from ase.utils import import_module
from ase.db.core import parallel, parallel_generator

all_formats = [
    'abinit', 'aims', 'aims_output', 'bundletrajectory', 'castep', 'castep_cell',
    'castep_geom', 'cfg', 'cif', 'cmdft', 'cube', 'dacapo', 'dacapo_text',
    'db', 'dftb', 'eon', 'eps', 'espresso_in', 'espresso_out', 'etsf', 'exciting', 'extxyz',
    'findsym', 'gromos', 'gaussian', 'gaussian_out', 'gen', 'gpaw_out', 'gpw',
    'gromacs', 'html', 'iwm', 'json', 'lammps_dump', 'mol', 'nwchem', 'pdb', 'png',
    'postgresql', 'pov', 'res', 'sdf', 'struct', 'struct_out', 'turbomole',
    'turbomole_gradient', 'traj', 'trj', 'v_sim', 'vasp', 'vasp_out',
    'vasp_xdatcar', 'vasp_xml', 'vti', 'vts', 'x3d', 'xsd', 'xsf', 'xyz']

if 0:
    import textwrap
    print(textwrap.fill("'" + "', '".join(sorted(set(all_formats))) + "'",
                        initial_indent='    ',
                        subsequent_indent='    ',
                        width=79))

# Special cases:
format2modulename = dict(
    traj='trajectory',
    json='db',
    postgresql='db',
    struct='wien2k',
    vti='vtkxml',
    vts='vtkxml',
    castep_cell='castep',
    castep_geom='castep',
    aims_out='aims',
    dacapo_text='dacapo',
    espresso_in='espresso',
    espresso_out='espresso',
    aims_output='aims',
    gaussian_out='gaussian',
    lammps_dump='lammpsrun',
    struct_out='siesta',
    turbomole_gradient='turbomole',
    trj='pickletrajectory',
    vasp_out='vasp',
    vasp_xdatcar='vasp',
    vasp_xml='vasp',
    html='x3d')

extension2format = dict(
    shelx='res',
    con='eon',
    cell='castep_cell',
    geom='castep_geom',
    out='espresso_out',
    gro='gromacs',
    g96='gromos',
    log='gaussian_out',
    com='gaussian',
    nw='nwchem',
    exi='exciting')

stores_multiple_images = [
    'xyz', 'traj', 'trj', 'pdb', 'cif', 'extxyz', 'db', 'json',
    'postgresql', 'xsf', 'findsym', 'gpaw_out', 'turbomole_gradient']

does_not_accept_a_file_descriptor = [
    'traj', 'db', 'postgresql',
    'etsf', 'dftb', 'aims', 'bundletrajectory', 'castep_cell', 'struct', 'res', 'eps', 'gpaw_out', 'gromacs', 'x3d', 'pov', 'trj', 'html']

IOFormat = collections.namedtuple('IOFormat', 'read, write, single, acceptsfd')
ioformats = {}  # will be filled at run-time
        
        
def initialize(format):
    if format in ioformats:
        return
    module_name = format2modulename.get(format, format)
    try:
        module = import_module('ase.io.' + module_name)
    except ImportError:
        print(format, module_name)
        raise ValueError('File format not recognized: ' + format)
    read = getattr(module, 'read_' + format, None)
    write = getattr(module, 'write_' + format, None)
    if read and not inspect.isgeneratorfunction(read):
        read = functools.partial(wrap_old_read_function, read)
    if not read and not write:
        raise ValueError('File format not recognized: ' + format)
    ioformats[format] = IOFormat(
        read, write,
        format not in stores_multiple_images,
        format not in does_not_accept_a_file_descriptor)
    

def get_ioformat(format):
    initialize(format)
    return ioformats[format]
    

def wrap_old_read_function(read, filename, index=None, **kwargs):
    if index is None:
        yield read(filename, **kwargs)
    else:
        for atoms in read(filename, index, **kwargs):
            yield atoms
        
        
@parallel
def write(filename, images, format=None, **kwargs):
    """Write Atoms object(s) to file.

    filename: str
        Name of the file to write to.
    images: Atoms object or list of Atoms objects
        A single Atoms object or a list of Atoms objects.
    format: str
        Used to specify the file-format.  If not given, the
        file-format will be taken from suffix of the filename.

    The accepted output formats:

    Many formats allow on open file-like object to be passed instead
    of ``filename``. In this case the format cannot be auto-decected,
    so the ``format`` argument should be explicitly given.

    The use of additional keywords is format specific."""

    if isinstance(filename, str):
        fd = None
        if filename == '-':
            fd = sys.stdout
            filename = None
        elif format is None:
            format = filetype(filename, read=False)
    else:
        fd = filename
        filename = None
        
    format = format or 'json'

    io = get_ioformat(format)

    if isinstance(images, Atoms):
        images = [images]
        
    if io.single:
        if len(images) > 1:
            raise ValueError('{0}-format can only store 1 Atoms object.'
                             .format(format))
        images = images[0]
        
    if io.write is None:
        raise ValueError("Can't write to {0}-format".format(format))
        
    # Special case for json-format:
    if format == 'json' and len(images) > 1:
        if filename is not None:
            io.write(filename, images, **kwargs)
            return
        raise ValueError("Can't write more than one image to file-descriptor"
                         'using json-format.')
        
    if io.acceptsfd:
        if fd is None:
            fd = open(filename, 'w')
        io.write(fd, images, **kwargs)
        fd.close()
    else:
        if fd is not None:
            raise ValueError("Can't write {0}-format to file-descriptor"
                             .format(format))
        io.write(filename, images, **kwargs)
    
    
def read(filename, index=None, format=None, **kwargs):
    """Read Atoms object(s) from file.

    filename: str
        Name of the file to read from.
    index: int or slice
        If the file contains several configurations, the last configuration
        will be returned by default.  Use index=n to get configuration
        number n (counting from zero).
    format: str
        Used to specify the file-format.  If not given, the
        file-format will be guessed by the *filetype* function.
        
    Many formats allow on open file-like object to be passed instead
    of ``filename``. In this case the format cannot be auto-decected,
    so the ``format`` argument should be explicitly given."""

    if isinstance(index, str):
        index = string2index(index)
    filename, index = parse_filename(filename, index)
    if index is None:
        index = -1
    if isinstance(index, (slice, str)):
        return list(_iread(filename, index, format, **kwargs))
    else:
        return next(_iread(filename, slice(index, None), format, **kwargs))
    
        
def iread(filename, index=None, format=None, **kwargs):
    if isinstance(index, str):
        index = string2index(index)
        
    filename, index = parse_filename(filename, index)
    
    if index is None:
        index = slice(None, None, None)
        
    if not isinstance(index, (slice, str)):
        index = slice(index, (index + 1) or None)
        
    for atoms in _iread(filename, index, format, **kwargs):
        yield atoms

            
@parallel_generator
def _iread(filename, index, format, **kwargs):
    if format is None:
        format = filetype(filename)

    io = get_ioformat(format)
    
    if not io.read:
        raise ValueError("Can't read from {0}-format".format(format))
        
    if io.single:
        start = index.start
        assert start is None or start == 0 or start == -1
        args = ()
    else:
        args = (index,)
        
    if isinstance(filename, str):
        if io.acceptsfd:
            fd = open(filename)
        else:
            fd = filename
    else:
        assert io.acceptsfd
        fd = filename
        
    # Make sure fd is closed in case loop doesn't finish:
    try:
        for atoms in io.read(fd, *args, **kwargs):
            yield atoms
    finally:
        if not isinstance(fd, str):
            fd.close()
    
    
def parse_filename(filename, index):
    if not isinstance(filename, str) or '@' not in filename:
        return filename, index
    newindex = None
    if ('.json@' in filename or
        '.db@' in filename or
        filename.startswith('pg://')):
        newfilename, newindex = filename.rsplit('@', 1)
    else:
        newfilename, newindex = filename.rsplit('@', 1)
        try:
            newindex = string2index(newindex)
        except ValueError:
            return filename, index
    if index is not None:
        raise ValueError('Only one index is allowed')
    return newfilename, newindex


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


def filetype(filename, read=True):
    """Try to guess the type of the file."""
    if isinstance(filename, str):
        if os.path.isdir(filename):
            if os.path.basename(os.path.normpath(filename)) == 'states':
                return 'eon'
            return 'bundletrajectory'

        if filename.startswith('pg://'):
            return 'postgresql'
        
        basename = os.path.basename(filename)
        
        if '.' in basename:
            ext = filename.rsplit('.', 1)[-1].lower()
            if ext in ['xyz', 'cube', 'json']:
                return ext

        if 'POSCAR' in basename or 'CONTCAR' in basename:
            return 'vasp'
        if 'OUTCAR' in basename:
            return 'vasp_out'
        if 'XDATCAR' in basename:
            return 'vasp_xdatcar'
        if 'vasp' in basename and basename.endswith('.xml'):
            return 'vasp_xml'
        if basename == 'coord':
            return 'turbomole'
        if basename == 'gradient':
            return 'turbomole_gradient'
        if basename.endswith('I_info'):
            return 'cmdft'
        if basename == 'atoms.dat':
            return 'iwm'
            
        if not read:
            return extension2format.get(ext, ext)
    
        fd = open(filename, 'rb')
    else:
        ext = 'json'
        fd = filename
        if fd is sys.stdin:
            return 'json'
            
    data = fd.read(2000)
    if fd is not filename:
        fd.close()
        
    if len(data) == 0:
        raise IOError('Empty file: ' + filename)

    for format, magic in [('traj', b'AFFormatASE-Trajectory'),
                          ('trj', b'PickleTrajectory'),
                          ('etsf', b'CDF'),
                          ('turbomole', b'$coord'),
                          ('turbomole_gradient', b'$grad'),
                          ('dftb', b'Geometry')]:
        if data.startswith(magic):
            return format

    for format, magic in [('gpaw-text', b'  ___ ___ ___ _ _ _  \n'),
                          ('esp_in', b'\n&system'),
                          ('esp_in', b'\n&SYSTEM'),
                          ('aims_out', b'\nInvoking FHI-aims ...'),
                          ('lammps_dump', b'\nITEM: TIMESTEP\n'),
                          ('xsf', b'\nANIMSTEPS'),
                          ('xsf', b'\nCRYSTAL'),
                          ('xsf', b'\nSLAB'),
                          ('xsf', b'\nPOLYMER'),
                          ('xsf', b'\nMOLECULE'),
                          ('xsf', b'\nATOMS'),
                          ('dacapo-text',
                           b'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n')]:
        if magic in data:
            return format

    return extension2format.get(ext, ext)
        
    
if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser()
    opts, filenames = parser.parse_args()
    n = max(len(filename) for filename in filenames) + 2
    for filename in filenames:
        try:
            format = filetype(filename)
        except ValueError:
            format = '?'
            
        # get description ...
        print('{0:{1}}{2}'.format(filename + ':', n, format))
