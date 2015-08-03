import collections
import contextlib
import functools
import inspect
import os
import sys

from ase.atoms import Atoms
from ase.utils import import_module
from ase.db.core import parallel, parallel_generator

format2modulename = dict(
    traj='trajectory',
    json='db',
    postgresql='db',
    cell='castep',
    exi='exciting',
    tmol='turbomole',
    struct='wien2k',
    gro='gromacs',
    g96='gromos',
    html='x3d')

extension2format = dict(
    shelx='res',
    con='eon')

not_single = ['xyz', 'traj', 'trj', 'pdb', 'cif', 'extxyz', 'db', 'json',
              'postgresql', 'xsf', 'findsym']
not_acceptsfd = ['traj', 'db', 'postgresql',
                 'struct', 'res', 'eps']

IOFormat = collections.namedtuple('IOFormat', 'read, write, single, acceptsfd')
# Will be filled at run-time:
ioformats = {}
        
        
def initialize(format):
    if format in ioformats:
        return
    module_name = format2modulename.get(format, format)
    try:
        module = import_module('ase.io.' + module_name)
    except ImportError:
        raise ValueError('File format not recognized: ' + format)
    read = getattr(module, 'read_' + format, None)
    write = getattr(module, 'write_' + format, None)
    if read and not inspect.isgeneratorfunction(read):
        read = functools.partial(convert_old_read_function, read)
    if not read and not write:
        raise ValueError('File format not recognized: ' + format)
    ioformats[format] = IOFormat(read, write,
                                 format not in not_single,
                                 format not in not_acceptsfd)
    

def get_ioformat(format):
    initialize(format)
    return ioformats[format]
    

def convert_old_read_function(read, filename, index=None, **kwargs):
    if index is None:
        yield read(filename, **kwargs)
    else:
        images = read(filename, index, **kwargs)
        #if isinstance(images, Atoms):
            #images = [images]
        for atoms in images:
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

    fd = None
    if isinstance(filename, str):
        if filename == '-':
            fd = sys.stdout
        elif format is None:
            format = filetype(filename, read=False)
    else:
        fd = filename
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
    if not io.acceptsfd:
        if fd is not None:
            raise ValueError("Can't write {0}-format to file-descriptor"
                             .format(format))
        io.write(filename, images, **kwargs)
    else:
        if fd is None:
            fd = open(filename, 'w')
        io.write(fd, images, **kwargs)
        fd.close()
    
    
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
        return list(iread(filename, index, format, **kwargs))
    else:
        return next(iread(filename, slice(index, None), format, **kwargs))
    
        
@parallel_generator
def iread(filename, index=None, format=None, **kwargs):
    if isinstance(index, str):
        index = string2index(index)
        
    filename, index = parse_filename(filename, index)
    
    if index is None:
        index = slice()
        
    if not isinstance(index, (slice, str)):
        index = slice(index, (index + 1) or None)
        
    if format is None:
        format = filetype(filename)

    io = get_ioformat(format)
    
    if not io.read:
        raise ValueError("Can't read from {0}-format".format(format))
        
    if io.single:
        start = index.start
        assert start is None or start == 0 or start == -1
        index = None
        
    if isinstance(filename, str):
        if io.acceptsfd:
            fd = open(filename)
        else:
            for atoms in io.read(filename, index, **kwargs):
                yield atoms
            return
    else:
        assert io.acceptsfd
        fd = filename
    with contextlib.closing(fd):
        for atoms in io.read(fd, index, **kwargs):
            yield atoms
    
    
def parse_filename(filename, index):
    if not isinstance(filename, str) or '@' not in filename:
        return filename, index
    newindex = None
    if ('.json@' in filename or '.db@' in filename or
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
            from ase.io.bundle import BundleTrajectory
            if BundleTrajectory.is_bundle(filename):
                return 'bundle'
            elif os.path.basename(os.path.normpath(filename)) == 'states':
                return 'eon'
            else:
                raise IOError('Directory: ' + filename)

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
            return 'tmol'
        if basename == 'gradient':
            return 'tmol-gradient'

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
                          ('coord', b'$coord'),
                          ('tmol', b'$grad'),
                          ('dftb', b'Geometry')]:
        if data.startswith(magic):
            return format

    for format, magic in [('gpaw-text', b'  ___ ___ ___ _ _ _  \n'),
                          ('esp_in', b'\n&system'),
                          ('esp_in', b'\n&SYSTEM'),
                          ('aims_out', b'\nInvoking FHI-aims ...'),
                          ('lammps', b'\nITEM: TIMESTEP\n'),
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

    return ext
        
    
if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser()
    opts, filenames = parser.parse()
    n = max(len(filename) for filename in filenames) + 2
    for filename in filenames:
        try:
            format = filetype(filename)
        except ValueError:
            format = '?'
            
        # get description ...
        print('{0:{1}}{2}'.format(filename + ':', n, format))
