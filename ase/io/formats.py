import contextlib
import inspect
import os
import sys

from ase.utils import import_module
from ase.atoms import Atoms
# from ase.parallel import run_generator_on_master, run_function_on_master

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

# Will be filled at run-time:
format_dicts = {}
        
        
def initialize(format):
    if format in format_dicts:
        return
    module_name = format2modulename.get(format, format)
    try:
        module = import_module('ase.io.' + module_name)
    except ImportError:
        raise ValueError('File format not recognized: ' + format)
    read = getattr(module, 'read_' + format, None)
    write = getattr(module, 'write_' + format, None)
    dct = {}
    if write:
        dct['write'] = write
    if read:
        dct['read'] = read
    if not dct:
        raise ValueError('File format not recognized: ' + format)
    format_dicts[format] = dct
    
        
#@run_function_on_master
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
    if isinstance(images, Atoms):
        images = [images]
    fd = None
    if isinstance(filename, str):
        if filename == '-':
            fd = sys.stdout
        elif format is None:
            format = filetype(filename, read=False)
    else:
        fd = filename
    format = format or 'json'
    initialize(format)
    if len(images) > 1 and format not in ['traj', 'xyz']:
        raise ValueError('{0}-format con only store 1 Atoms object.'
                         .format(format))
    write = format_dicts[format].get('write')
    if write is None:
        raise ValueError("Can't write to {0}-format".format(format))
    if format in ['traj', 'db']:
        if fd is not None:
            raise ValueError("Can't write {0}-format to file-descriptor"
                             .format(format))
        write(filename, images, **kwargs)
    else:
        if fd is None:
            fd = open(filename, 'w')
        write(fd, images, **kwargs)
        fd.close()
    
    
def read(filename, index=None, format=None):
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
        return list(iread(filename, index, format))
    else:
        return next(iread(filename, slice(index, None), format))
    
        
#@run_generator_on_master
def iread(filename, index=None, format=None):
    if isinstance(index, str):
        index = string2index(index)
        
    filename, index = parse_filename(filename, index)
    
    if isinstance(index, str):
        import ase.db
        con = ase.db.connect(filename)
        for row in con.select(index):
            yield row.toatoms()
        return
        
    if index is None:
        index = slice()
        
    if not isinstance(index, (slice, str)):
        index = slice(index, (index + 1) or None)
        
    if format is None:
        format = filetype(filename)
    initialize(format)
    read = format_dicts[format].get('read')
    
    if not read:
        raise ValueError("Can't read from {0}-format".format(format))
        
    if inspect.isgeneratorfunction(read):
        if isinstance(filename, str):
            fd = open(filename)
        else:
            fd = filename
            
        with contextlib.closing(fd):
            for atoms in read(fd, index):
                yield atoms
    else:
        images = read(filename, index)
        if isinstance(images, Atoms):
            images = [images]
        for atoms in images:
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
            if ext in ['xyz', 'cube']:
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
