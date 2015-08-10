# coding=utf-8
"""ASE-interface to Octopus.

Ask Hjorth Larsen <asklarsen@gmail.com>
Carlos de Armas

http://tddft.org/programs/octopus/
"""
import os

import numpy as np

from ase import Atoms
from ase.calculators.calculator import FileIOCalculator
# XXX raise ReadError upon bad read
from ase.data import atomic_numbers
from ase.io import read
from ase.io.xsf import read_xsf

# Representation of parameters from highest to lowest level of abstraction:
#
#  * Atoms object plus reduced kwargs that specify info not stored in the Atoms
#  * full dictionary of kwargs incorporating info contained in the Atoms
#  * series of assignments (names, values).  May contain duplicates.
#  * text in Octopus input file format


# Octopus variable types and specification in Python:
#
#  Type     Examples                    Equivalent in Python:
# -----------------------------------------------------------------------------
#  flag     wfs + density               'wfs + density'
#  float    2.7, 2.7 + pi^2             2.7, 2.7 + np.pi**2
#  integer  42, rmmdiis                 42, 'rmmdiis'
#  logical  true, false, yes, no, ...   True, False, 1, 0, 'yes', 'no', ...
#  string   "stdout.txt"                '"stdout.txt"' (apologies for ugliness)
#
#  block    %Coordinates                List of lists:
#            'H' | 0 | 0 | 0              coordinates=[["'H'", 0, 0, 0],
#            'O' | 0 | 0 | 1                           ["'O'", 0, 0, 1]]
#           %                             (elemements are sent through repr())

# Rules for input parameters
# --------------------------
#
# We make the following conversions:
#     dict of keyword arguments + Atoms object -> Octopus input file
# and:
#     Octopus input file -> Atoms object + dict of keyword arguments
# Below we detail some conventions and compatibility issues.
#
# 1) ASE always passes some parameters by default (Units=eV_Angstrom,
#    etc.).  They can be overridden by the user but the resulting
#    behaviour is undefined.
#
# 2) Atoms object is used to establish some parameters: Coordinates,
#    Lsize, etc.  All those parameters can be overridden by passing
#    them directly as keyword arguments.  Parameters that were taken
#    from the Atoms object are always marked with the comment "# ASE
#    auto" in the input file.  This is used to distinguish variables
#    that are overridden from variables that simply came from the
#    atoms object when restarting.
#
# 3) Some variables do not interact nicely between ASE and Octopus,
#    such as SubSystemCoordinates which may involve rotations.  There
#    may be many such variables that we have not identified, but at
#    least the known ones will cause a suppressable
#    OctopusKeywordError.  (This third rule has not been implemented
#    as of this moment.)

class OctopusKeywordError(ValueError):
    pass  # Unhandled keywords


class OctopusParseError(Exception):
    pass  # Cannot parse input file


class OctopusIOError(IOError):
    pass  # Cannot find output files


def unpad(pbc, arr):
    # Return non-padded array from padded array.
    # This means removing the last element along all periodic directions.
    if pbc[0]:
        assert np.all(arr[0, :, :] == arr[-1, :, :])
        arr = arr[0:-1, :, :]
    if pbc[1]:
        assert np.all(arr[:, 0, :] == arr[:, -1, :])
        arr = arr[:, 0:-1, :]
    if pbc[2]:
        assert np.all(arr[:, :, 0] == arr[:, :, -1])
        arr = arr[:, :, 0:-1]
    return np.ascontiguousarray(arr)

    
def unpad_smarter(pbc, arr):
    # 'Smarter' but less easy to understand version of the above.
    # (untested I think)
    slices = []
    for c, is_periodic in enumerate(pbc):
        if is_periodic:
            left = np.take(arr, [0], axis=c)
            right = np.take(arr, [-1], axis=c)
            assert np.all(left == right)
            slices.append(slice(0, -1))
        else:
            slices.append(slice(None))
    return np.ascontiguousarray(arr[slices])


# Octopus writes slightly broken XSF files.  This hack purports to fix them.
def repair_brokenness_of_octopus_xsf(path):
    assert os.path.isfile(path), path

    def replace(old, new):
        # XXX Ouch... need to get rid of this
        os.system('sed -i s/%s/%s/ %s' % (old, new, path))
        
    replace('BEGIN_BLOCK_DATAGRID3D', 'BEGIN_BLOCK_DATAGRID_3D')
    replace('^DATAGRID_3D', 'BEGIN_DATAGRID_3D')


def fix_and_read_xsf(fname, read_data=False):
    repair_brokenness_of_octopus_xsf(fname)
    return read_xsf(fname, read_data=read_data)


# Parse value as written in input file *or* something that one would be
# passing to the ASE interface, i.e., this might already be a boolean
def octbool2bool(value):
    value = value.lower()
    if isinstance(value, int):
        return bool(value)
    if value in ['true', 't', 'yes', '1']:
        return True
    elif value in ['no', 'f', 'false', '0']:
        return False
    else:
        raise ValueError('Failed to interpret "%s" as a boolean.' % value)


def list2block(name, rows):
    """Construct 'block' of Octopus input.

    convert a list of rows to a string with the format x | x | ....
    for the octopus input file
    :rtype : string
    :param name: name of the block
    :param rows: list of data to put in the string
    :return: the format string"""
    lines = []
    lines.append('%' + name)
    for row in rows:
        lines.append(' ' + ' | '.join(str(obj) for obj in row))
    lines.append('%')
    return lines


def purify(kwargs):
    """Reduce keywords to unambiguous form (lowercase).

    Also convert to string all numbers.
    :param kwargs: parameter pass to Calculator class"""
    newkwargs = {}
    for arg, value in kwargs.items():
        lkey = arg.lower()
        newkwargs[lkey] = value
        if isinstance(value, int):
            newkwargs[lkey] = str(value)
        elif isinstance(value, float):
            newkwargs[lkey] = str(value)
    return newkwargs


def input_line_iter(lines):
    """Convenient iterator for parsing input files 'cleanly'.

    Discards comments etc."""
    for line in lines:
        line = line.split('#')[0].strip()
        if not line or line.isspace():
            continue
        line = line.strip()
        yield line


def block2list(lines, header=None):
    """Parse lines of block and return list of lists of strings."""
    lines = iter(lines)
    block = []
    if header is None:
        header = lines.next()
    assert header.startswith('%'), header
    name = header[1:]
    for line in lines:
        if line.startswith('%'):  # Could also say line == '%' most likely.
            break
        tokens = [token.strip() for token in line.strip().split('|')]
        block.append(tokens)
    return name, block


def parse_input_file(fd):
    names = []
    values = []
    lines = input_line_iter(fd)
    while True:
        try:
            line = lines.next()
        except StopIteration:
            break
        else:
            if line.startswith('%'):
                name, value = block2list(lines, header=line)
            else:
                tokens = line.split('=')
                assert len(tokens) == 2, tokens
                name = tokens[0].strip()
                value = tokens[1].strip()
            names.append(name)
            values.append(value)
    return names, values


def kwargs2cell(kwargs):
    # kwargs -> cell + remaining kwargs
    # cell will be None if not ASE-compatible.
    kwargs = purify(kwargs)

    if boxshape_is_ase_compatible(kwargs):
        kwargs.pop('boxshape', None)
        Lsize = kwargs.pop('lsize')
        assert len(Lsize) == 1
        cell = [2 * float(l) for l in Lsize[0]]
        # TODO support LatticeVectors
    else:
        cell = None
    return cell, kwargs


def boxshape_is_ase_compatible(kwargs):
    boxshape = kwargs.get('boxshape')
    # XXX add support for experimental keyword 'latticevectors'
    if boxshape == 'parallelepiped':
        return True
    else:
        pdims = kwargs.get('periodicdimensions', 0)
        return boxshape is None and pdims > 0


def kwargs2atoms(kwargs):
    """Extract atoms object from keywords and return remaining keywords."""
    kwargs = purify(kwargs)

    coord_keywords = ['coordinates',
                      'xyzcoordinates',
                      'pdbcoordinates',
                      'reducedcoordinates',
                      'xsfcoordinates',
                      'xsfcoordinatesanimstep']
    
    nkeywords = 0
    for keyword in coord_keywords:
        if keyword in kwargs:
            nkeywords += 1
    if nkeywords == 0:
        raise OctopusParseError('No coordinates')
    elif nkeywords > 1:
        raise OctopusParseError('Multiple coordinate specifications present.  '
                                'This may be okay in Octopus, but we do not '
                                'implement it.')

    def get_positions_from_block(keyword):
        # %Coordinates or %ReducedCoordinates -> atomic numbers, positions.
        block = kwargs.pop(keyword)
        positions = []
        numbers = []
        for row in block:
            assert len(row) in [ndims + 1, ndims + 2]
            row = row[:ndims + 1]
            sym = row[0]
            assert sym.startswith("'") and sym.endswith("'")
            sym = sym[1:-1]
            pos0 = np.zeros(3)
            ndim = kwargs.get('dimensions', 3)
            pos0[:ndim] = [float(element) for element in row[1:]]
            number = atomic_numbers[sym]  # Use 0 ~ 'X' for unknown?
            numbers.append(number)
            positions.append(pos0)
        positions = np.array(positions)
        return numbers, positions

    def read_atoms_from_file(keyword):
        if keyword not in kwargs:
            return None
        
        fname = kwargs.pop(keyword)
        fmt = keyword[:3].lower()
        # XXX test xyz, pbd and xsf
        if fmt == 'xsf' and 'xsfcoordinatesanimstep' in kwargs:
            anim_step = kwargs.pop('xsfcoordinatesanimstep')
            slice = slice(anim_step, anim_step + 1, 1)
            # XXX test animstep
        images = read(fname, slice(None, None, 1), fmt)
        if len(images) != 1:
            raise OctopusParseError('Expected only one image.  Don\'t know '
                                    'what to do with %d images.' % len(images))
        return images[0]

    ndims = kwargs.get('dimensions', 3)
    cell, kwargs = kwargs2cell(kwargs)
    # XXX fix interaction between pbc and possibly existing pbc in XSF/etc.

    # XXX remember to pop
    atoms = None
    if 'coordinates' in kwargs:
        numbers, positions = get_positions_from_block('coordinates')
        if cell is not None:
            positions[:, :] += np.array(cell)[None, :] / 2.0
        atoms = Atoms(cell=cell, numbers=numbers, positions=positions)
    elif 'reducedcoordinates' in kwargs:
        numbers, rpositions = get_positions_from_block('reducedcoordinates')
        if cell is None:
            raise ValueError('Cannot figure out what the cell is, '
                             'and thus cannot interpret reduced coordinates.')
        rpositions += 0.5
        atoms = Atoms(cell=cell, numbers=numbers, scaled_positions=rpositions)
    else:
        for keyword in ['xyzcoordinates', 'pdbcoordinates', 'xsfcoordinates']:
            atoms = read_atoms_from_file(keyword)
            if atoms is not None:
                break
        else:
            raise OctopusParseError('Apparently there are no atoms.')

    assert atoms is not None

    # Either we have non-periodic BCs or the atoms object already
    # got its BCs from reading the file.  In the latter case
    # we shall override only if PeriodicDimensions was given specifically:
    pdims = kwargs.pop('periodicdimensions', None)
    if pdims is not None:
        pbc = np.zeros(3, dtype=bool)
        pdims = int(pdims)
        pbc[:pdims] = True
        atoms.pbc = pbc

    if cell is not None:
        atoms.cell = cell

    return atoms, kwargs


def atoms2kwargs(atoms):
    kwargs = {}

    kwargs['boxshape'] = 'parallelepiped'
    # TODO LatticeVectors parameter for non-orthogonal cells
    Lsize = 0.5 * np.diag(atoms.cell).copy()
    kwargs['lsize'] = [[repr(size) for size in Lsize]]

    # ASE uses (0...cell) while Octopus uses -L/2...L/2.
    # Lsize is really cell / 2, and we have to adjust our
    # positions by subtracting Lsize (see construction of the coords
    # block).
    positions = atoms.positions.copy()
    positions -= Lsize[None, :]

    coord_block = [[repr(sym)] + list(map(repr, pos))
                   for sym, pos in zip(atoms.get_chemical_symbols(),
                                       positions)]
    kwargs['coordinates'] = coord_block
    npbc = sum(atoms.pbc)
    for c in range(npbc):
        if not atoms.pbc[c]:
            msg = ('Boundary conditions of Atoms object inconsistent '
                   'with requirements of Octopus.  pbc must be either '
                   '000, 100, 110, or 111.')
            raise ValueError(msg)
    kwargs['periodicdimensions'] = npbc

    # TODO InitialSpins
    #
    # TODO can use maximumiterations + output/outputhow to extract
    # things from restart file into output files without trouble.
    #
    # Velocities etc.?
    return kwargs


def generate_input(atoms, kwargs):
    """Convert atoms and keyword arguments to Octopus input file."""
    _lines = []

    def append(line):
        _lines.append(line)

    def extend(lines):
        _lines.extend(lines)
        append('')

    def setvar(key, var):
        append('%s = %s' % (key, var))

    if 'units' in kwargs:
        raise ValueError('Sorry, but we decide the units in the ASE '
                         'interface for now.')
    setvar('units', 'eV_Angstrom')
    
    atomskwargs = atoms2kwargs(atoms)

    # Use cell from Atoms object unless user specified BoxShape
    use_ase_box = 'boxshape' not in kwargs
    if use_ase_box:
        setvar('boxshape', atomskwargs['boxshape'])
        lsizeblock = list2block('lsize', atomskwargs['lsize'])
        extend(lsizeblock)

    # Allow override or issue errors?
    pdim = 'periodicdimensions'
    if pdim in kwargs:
        if int(kwargs[pdim]) != int(atomskwargs[pdim]):
            raise ValueError('Cannot reconcile periodicity in input '
                             'with that of Atoms object')
    setvar('periodicdimensions', atomskwargs[pdim])

    # We like to output forces
    if 'output' in kwargs:
        output_string = kwargs.pop('output')
        output_tokens = [token.strip()
                         for token in output_string.lower().split('+')]
    else:
        output_tokens = []

    if 'forces' not in output_tokens:
        output_tokens.append('forces')
    setvar('output', ' + '.join(output_tokens))
    # It is illegal to have output forces without any OutputHow.
    # Even though the forces are written in the same format no matter
    # OutputHow.  Thus we have to make one up:
    if 'outputhow' not in kwargs:
        kwargs['outputhow'] = 'xcrysden'

    for key in kwargs:
        val = kwargs[key]
        # Most datatypes are straightforward but blocks require some attention.
        if isinstance(val, list):
            append('')
            dict_data = list2block(key, val)
            extend(dict_data)
        else:
            setvar(key, str(val))
    append('')

    coord_block = list2block('Coordinates', atomskwargs['coordinates'])
    extend(coord_block)
    
    return '\n'.join(_lines)


class Octopus(FileIOCalculator):
    """Octopus calculator.

    The label is always assumed to be a directory."""

    implemented_properties = ['energy', 'forces',
                              'dipole', #'charges',
                              'magmom', 'magmoms']
    # valores propios en las opciones output y outputhow

    troublesome_keywords = set(['subsystemcoordinates',
                                'subsystems',
                                'unitsinput',
                                'unitsoutput',
                                'pdbcoordinates',
                                'xyzcoordinates',
                                'xsfcoordinates',
                                'xsfcoordinatesanimstep',
                                'reducedcoordinates'])

    def __init__(self,
                 restart=None,
                 label=None,
                 atoms=None,
                 command='octopus',
                 ignore_troublesome_keywords=None,
                 **kwargs):
        """Create Octopus calculator.

        Label is always taken as a subdirectory.
        Restart is taken to be a label."""

        # XXX support the specially defined ASE parameters,
        # "smear" etc.

        if restart is not None:
            if label is not None and restart != label:
                raise ValueError('restart and label are mutually exclusive '
                                 'or must at the very least coincide.')
            label = restart

        if label is None:
            label = 'ink-pool'

        if ignore_troublesome_keywords:
            trouble = set(self.troublesome_keywords)
            for keyword in ignore_troublesome_keywords:
                trouble.pop(keyword)
            # XXX test this
            self.troublesome_keywords = trouble

        self.kwargs = {}

        FileIOCalculator.__init__(self, restart=restart,
                                  ignore_bad_restart_file=False,
                                  label=label,
                                  atoms=atoms,
                                  command=command, **kwargs)
        # The above call triggers set() so we can update self.kwargs.

    def set_label(self, label):
        # Octopus does not support arbitrary namings of all the output files.
        # But we can decide that we always dump everything in a directory.
        if not label.endswith('/'):
            label += '/'
        FileIOCalculator.set_label(self, label)

    def set(self, **kwargs):
        """Set octopus input file parameters.
        
        :param kwargs:
        """
        kwargs = purify(kwargs)
        for keyword in kwargs:
            if keyword in self.troublesome_keywords:
                msg = ('ASE-Octopus interface will probably misbehave with '
                       'the %s parameter.  Optimists may use '
                       'Octopus(ignore_troublesome_keywords=[kw1, kw2, ...])'
                       'to override this.' % keyword)
                raise OctopusKeywordError(msg)

        FileIOCalculator.set(self, **kwargs)
        self.kwargs.update(kwargs)
        # XXX should use 'Parameters' but don't know how

    def get_xc_functional(self):
        """Return the XC-functional identifier.
            'LDA', 'PBE', ..."""
        return self.kwargs.get('xcfunctional', 'LDA')

    def get_bz_k_points(self):
        """Return all the k-points in the 1. Brillouin zone.
        The coordinates are relative to reciprocal latice vectors."""
        # Have not found nice way of extracting this information
        # from Octopus.  Thus unimplemented. -askhl
        raise NotImplementedError

    def get_charges(self, atoms=None):
        raise NotImplementedError

    def get_fermi_level(self):
        return self.results['efermi']

    def get_potential_energies(self):
        raise NotImplementedError

    def get_dipole_moment(self, atoms=None):
        if not 'dipole' in self.results:
            msg = ('Dipole moment not calculated.\n'
                   'You may wish to use SCFCalculateDipole=True')
            raise OctopusIOError(msg)
        return self.results['dipole']

    def get_stresses(self):
        raise NotImplementedError

    def _read_array(self, fname, outputkeyword=None):
        path = self._getpath('static/%s' % fname)
        if not os.path.exists(path):
            msg = 'Path not found: %s' % path
            if outputkeyword is not None:
                msg += ('\nIt appears that the %s has not been saved.\n'
                        'Be sure to specify Output=\'%s\' in the input.'
                        % (outputkeyword, outputkeyword))
            raise OctopusIOError(msg)
        # If this causes an error now that the file exists, things are
        # messed up.  Then it is better that the error propagates as normal
        return fix_and_read_xsf(path, read_data=True)

    def read_vn(self, basefname, keywordname):
        static_dir = self._getpath('static')
        assert os.path.isdir(static_dir)

        if self.get_spin_polarized():
            spin1, _atoms = self._read_array('%s-sp1.xsf' % basefname,
                                             keywordname)
            spin2, _atoms = self._read_array('%s-sp2.xsf' % basefname,
                                             keywordname)
            array = np.array([spin1, spin2])  # shape 2, nx, ny, nz
        else:
            array, _atoms = self._read_array('%s.xsf' % basefname, keywordname)
            array = array[None]  # shape 1, nx, ny, nx
        assert len(array.shape) == 4
        return array

    def _unpad_periodic(self, array):
        return unpad(self.get_atoms().pbc, array)

    def _pad_unperiodic(self, array):
        pbc = self.get_atoms().pbc
        orig_shape = array.shape
        newshape = [orig_shape[c] + (0 if pbc[c] else 1) for c in range(3)]
        out = np.zeros(newshape, dtype=array.dtype)
        nx, ny, nz = orig_shape
        out[:nx, :ny, :nz] = array
        return out

    def _pad_correctly(self, array, pad):
        array = self._unpad_periodic(array)
        if pad:
            array = self._pad_unperiodic(array)
        return array

    def get_pseudo_density(self, spin=None, pad=True):
        """Return pseudo-density array.

        If *spin* is not given, then the total density is returned.
        Otherwise, the spin up or down density is returned (spin=0 or
        1)."""
        if 'density_sg' not in self.results:
            self.results['density_sg'] = self.read_vn('density', 'density')
        density_sg = self.results['density_sg']
        if spin is None:
            density_g = density_sg.sum(axis=0)
        else:
            assert spin == 0 or (spin == 1 and len(density_sg) == 2)
            density_g = density_sg[spin]
        return self._pad_correctly(density_g, pad)

    def get_effective_potential(self, spin=0, pad=True):
        if spin is None:  # Annoying case because it works as an index!
            raise ValueError('spin=None')
        if 'potential_sg' not in self.results:
            self.results['potential_sg'] = self.read_vn('vks', 'potential')
        array = self.results['potential_sg'][spin]
        return self._pad_correctly(array, pad)

    def get_pseudo_wave_function(self, band=0, kpt=0, spin=0, broadcast=True,
                                 pad=True):
        """Return pseudo-wave-function array."""
        assert band < self.get_number_of_bands()

        ibz_k_pts = self.get_ibz_k_points()

        forcecomplex = self.kwargs.get('forcecomplex')
        if forcecomplex is not None:
            forcecomplex = octbool2bool(forcecomplex)
        if len(ibz_k_pts) > 1 or ibz_k_pts.any() or forcecomplex:
            dtype = complex
        else:
            dtype = float  # Might there be more issues that determine dtype?

        if self.get_spin_polarized():
            kpt_index = 2 * kpt + spin  # XXX this is *probably* correct
        else:
            kpt_index = kpt

        # The ASE convention is that kpts and bands start from 0,
        # whereas in Octopus they start from 1.  So always add 1
        # when looking for filenames.
        kpt_index += 1
        band_index = band + 1

        tokens = ['wf']
        if len(ibz_k_pts) > 1 or self.get_spin_polarized():
            tokens.append('-k%03d' % kpt_index)
        tokens.append('-st%04d' % band_index)
        name = ''.join(tokens)

        if dtype == float:
            array, _atoms = self._read_array('%s.xsf' % name, 'wfs')
        else:
            array_real, _atoms = self._read_array('%s.real.xsf' % name, 'wfs')
            array_imag, _atoms  = self._read_array('%s.imag.xsf' % name, 'wfs')
            array = array_real + 1j * array_imag

        return self._pad_correctly(array, pad)

    def get_number_of_spins(self):
        """Return the number of spins in the calculation.
           Spin-paired calculations: 1, spin-polarized calculation: 2."""
        return 2 if self.get_spin_polarized() else 1

    def get_spin_polarized(self):
        """Is it a spin-polarized calculation?"""

        sc = self.kwargs.get('spincomponents')
        if sc is None or sc == 'unpolarized':
            return False
        elif sc == 'spin_polarized' or sc == 'polarized':
            return True
        else:
            raise NotImplementedError('SpinComponents keyword %s' % sc)

    def get_ibz_k_points(self):
        """Return k-points in the irreducible part of the Brillouin zone.
        The coordinates are relative to reciprocal latice vectors."""
        return self.results['ibz_k_points']

    def get_k_point_weights(self):
        return self.results['k_point_weights']

    def get_number_of_bands(self):
        return self.results['nbands']

    def get_magnetic_moments(self, atoms=None):
        return self.results['magmoms'].copy()

    def get_magnetic_moment(self, atoms=None):
        return self.results['magmom']

    def get_occupation_numbers(self, kpt=0, spin=0):
        return self.results['occupations'][spin, kpt].copy()

    def get_eigenvalues(self, kpt=0, spin=0):
        return self.results['eigenvalues'][spin, kpt].copy()

    def _getpath(self, path, check=False):
        path = os.path.join(self.directory, path)
        if check:
            if not os.path.exists(path):
                raise OctopusIOError('No such file or directory: %s' % path)
        return path

    def get_atoms(self):
        return FileIOCalculator.get_atoms(self)

    def read_results(self):
        """Read octopus output files and extract data."""
        fd = open(self._getpath('static/info', check=True))

        def search(token):
            initial_pos = fd.tell()
            for line in fd:
                if line.strip().startswith(token):
                    return line
            try:
                raise ValueError('No such token: \'%s\'' % token)
            finally:
                fd.seek(initial_pos)

        try:
            nkptsline = search('Number of symmetry-reduced')
        except ValueError:
            nkpts = 1
            ibz_k_points = np.zeros((1, 3))
            k_point_weights = np.ones(1)
        else:
            nkpts = int(nkptsline.split('=')[-1].strip())

            kpt_coord_weight_line = search('List of k-points:')
            kpt_coord_weight_line = fd.next()
            tokens = kpt_coord_weight_line.split()
            assert tokens == ['ik', 'k_x', 'k_y', 'k_z', 'Weight']
            kpt_coord_weight_line = fd.next()
            kpt_coord_weight_line = fd.next()
            ibz_k_points = np.zeros(shape=[nkpts, 3])
            k_point_weights = np.zeros(shape=[nkpts])
            for i in range(nkpts):
                tokens = kpt_coord_weight_line.split()
                assert int(tokens[0])
                ibz_k_points[i][0] = float(tokens[1])
                ibz_k_points[i][1] = float(tokens[2])
                ibz_k_points[i][2] = float(tokens[3])
                k_point_weights[i] = float(tokens[4])
                kpt_coord_weight_line = fd.next()

        self.results['ibz_k_points'] = ibz_k_points
        self.results['k_point_weights'] = k_point_weights

        statesheader = search('Eigenvalues [')
        assert statesheader == 'Eigenvalues [eV]\n'

        for line in fd:
            tokens = line.split()
            try:
                if int(tokens[0]):
                    break
            except ValueError:
                continue

        kpts = []
        kpts_down = []

        class State(object):
            """Represents pairs (energy, occupation)."""
            def __init__(self, energy, occ):
                self.energy = energy
                self.occ = occ

        for k in range(nkpts):
            nbands = None
            states = []
            states_down = []

            if line.startswith('#'):
                line = fd.next()
            while line.strip(' ')[0].isdigit():
                tokens = line.split()
                # n = int(tokens[0])
                if tokens[1] in ['--', 'up']:
                    spin = 0
                elif tokens[1] == 'dn':
                    spin = 1
                else:
                    raise ValueError('Unexpected spin specification: %s'
                                     % tokens[1])
                energy = float(tokens[2])
                occupation = float(tokens[3])
                state = State(energy, occupation)
                if spin == 0:
                    states.append(state)
                else:
                    states_down.append(state)
                line = fd.next()
            if nbands is None:
                nbands = len(states)
                assert nbands > 0
            else:
                # this assert is triggered in the second run if
                # nbands is not reinitialized
                assert nbands == len(states)
            kpts.append(states)
            if len(states_down) > 0:
                kpts_down.append(states_down)

        nspins = 1
        if len(kpts_down):
            nspins = 2
            assert len(kpts_down) == len(kpts)
        eps_skn = np.empty((nspins, len(kpts), len(states)))
        occ_skn = np.empty((nspins, len(kpts), len(states)))

        for k in range(len(kpts)):
            kpt = kpts[k]
            eps_skn[0, k, :] = [state.energy for state in kpt]
            occ_skn[0, k, :] = [state.occ for state in kpt]
            if nspins == 2:
                kpt = kpts_down[k]
                eps_skn[1, k, :] = [state.energy for state in kpt]
                occ_skn[1, k, :] = [state.occ for state in kpt]

        self.results['nbands'] = nbands
        self.results['eigenvalues'] = eps_skn
        self.results['occupations'] = occ_skn

        if line.startswith('Fermi'):
            tokens = line.split()
            assert tokens[-1] == 'eV'
            eFermi = float(tokens[-2])
        else:
            # Find HOMO level.  Note: This could be a very bad
            # implementation with fractional occupations if the Fermi
            # level was not found otherwise.
            all_energies = eps_skn.ravel()
            all_occupations = occ_skn.ravel()
            args = np.argsort(all_energies)
            for arg in args[::-1]:
                if all_occupations[arg] > 0.1:
                    break
            eFermi = all_energies[arg]
        self.results['efermi'] = eFermi

        search('Energy [eV]:')
        line = fd.next()
        assert line.strip().startswith('Total'), line
        self.results['energy'] = float(line.split('=')[-1].strip())
        line = fd.next()
        assert line.strip().startswith('Free'), line
        self.results['free_energy'] = float(line.split('=')[-1].strip())

        if nspins == 2:
            line = search('Total Magnetic Moment:')
            line = fd.next()
            values = line.split()
            self.results['magmom'] = float(values[-1])

            line = fd.next()
            assert line.startswith('Local Magnetic Moments')
            line = fd.next()
            assert line.split() == ['Ion', 'mz']
            # Reading  Local Magnetic Moments
            mag_moment = []
            for line in fd:
                if line == '\n':
                    break  # there is no more thing to search for
                line = line.replace('\n', ' ')
                values = line.split()
                mag_moment.append(float(values[-1]))

            self.results['magmoms'] = np.array(mag_moment, dtype=float)
        else:
            self.results['magmom'] = 0.0
            self.results['magmoms'] = np.zeros(len(self.get_atoms()))

        # Read dipole data
        try:
            line = search('Dipole:')
        except ValueError:
            pass
        else:
            dipole = np.zeros(shape=[3, 2])
            line = fd.next()
            for i in range(3):
                line = line.replace('<', ' ')
                line = line.replace('>', ' ')
                line = line.replace('=', ' ')
                line = line.replace('\n', ' ')
                values = line.split()
                dipole[i][0] = float(values[1])
                dipole[i][1] = float(values[2])
                line = fd.next()

            self.results['dipole'] = dipole

        forces_atoms = read_xsf(self._getpath('static/forces.xsf'))
        F_av = forces_atoms.get_forces()
        self.results['forces'] = F_av

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties=properties,
                                     system_changes=system_changes)
        txt = generate_input(atoms, self.kwargs)
        fd = open(self._getpath('inp'), 'w')
        fd.write(txt)
        fd.close()

    def read(self, label):
        # XXX label of restart file may not be the same as actual label!
        # This makes things rather tricky.  We first set the label to
        # that of the restart file and arbitrarily expect the remaining code
        # to rectify any consequent inconsistencies.
        self.set_label(label)

        FileIOCalculator.read(self, label)
        inp_path = self._getpath('inp')
        fd = open(inp_path)
        names, values = parse_input_file(fd)
        kwargs = purify(dict(zip(names, values)))

        self.atoms, kwargs = kwargs2atoms(kwargs)
        self.kwargs.update(kwargs)

        fd.close()
        self.read_results()

    @classmethod
    def recipe(cls, **kwargs):
        system = Atoms()
        calc = Octopus(CalculationMode='recipe', **kwargs)
        system.set_calculator(calc)
        try:
            system.get_potential_energy()
        except OctopusIOError:
            pass
        else:
            raise OctopusIOError('Expected recipe, but found '
                                 'useful physical output!')

def main():
    from ase.lattice import bulk
    from ase.calculators.interfacechecker import check_interface

    system = bulk('Si', 'diamond', orthorhombic=True)
    calc = Octopus(Spacing=0.275,
                   KPointsGrid=[[2, 2, 2]],
                   KPointsUseSymmetries=True,
                   Smearing=0.1,
                   SmearingFunction='fermi_dirac',
                   ExtraStates=2,
                   stdout='"stdout.log"',
                   stderr='"stderr.log"',
                   Output='density + potential + wfs',
                   OutputHow='xcrysden')
    system.set_calculator(calc)
    system.get_potential_energy()

    check_interface(calc)

if __name__ == '__main__':
    main()
