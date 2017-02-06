"""Reads quantum espresso files. Tested for output on PWSCF v.5.0.2, only
for typical output of input files made with ASE -- that is, ibrav=0."""

import warnings
import operator as op
from collections import OrderedDict

import numpy as np
from ase.atoms import Atoms, Atom
from ase.units import create_units
from ase.calculators.singlepoint import SinglePointCalculator
from ase.utils import basestring, chemical_symbols

# Quantum ESPRESSO uses CODATA 2006 internally
units = create_units('2006')


def read_espresso_out(fileobj, index=-1):
    """Reads quantum espresso output text files."""
    if isinstance(fileobj, basestring):
        fileobj = open(fileobj, 'rU')
    lines = fileobj.readlines()
    images = []

    # Get unit cell info.
    bl_line = [line for line in lines if 'bravais-lattice index' in line]
    if len(bl_line) != 1:
        raise NotImplementedError('Unsupported: unit cell changing.')
    bl_line = bl_line[0].strip()
    brav_latt_index = bl_line.split('=')[1].strip()
    if brav_latt_index != '0':
        raise NotImplementedError('Supported only for Bravais-lattice '
                                  'index of 0 (free).')
    # get alat=celldm(1); use celldm it has more decimal places.
    lp_line = [line for line in lines if 'celldm(1)' in line]
    # TODO: implement changing cell shape
    if len(lp_line) != 1:
        raise NotImplementedError('Unsupported: unit cell changing.')
    lp_line = lp_line[0].split()[1]
    lattice_parameter = float(lp_line) * units['Bohr']
    ca_line_no = [number for (number, line) in enumerate(lines) if
                  'crystal axes: (cart. coord. in units of alat)' in line]
    if len(ca_line_no) != 1:
        raise NotImplementedError('Unsupported: unit cell changing.')
    ca_line_no = int(ca_line_no[0])
    cell = np.zeros((3, 3))
    for number, line in enumerate(lines[ca_line_no + 1: ca_line_no + 4]):
        line = line.split('=')[1].strip()[1:-1]
        values = [float(value) for value in line.split()]
        cell[number, 0] = values[0]
        cell[number, 1] = values[1]
        cell[number, 2] = values[2]
    cell *= lattice_parameter

    # Find atomic positions and add to images.
    for number, line in enumerate(lines):
        key = 'Begin final coordinates'  # these just reprint last posn.
        if key in line:
            break
        key = 'Cartesian axes'
        if key in line:
            atoms = make_atoms(number, lines, key, cell, lattice_parameter)
            images.append(atoms)
        key = 'ATOMIC_POSITIONS'
        if key in line:
            atoms = make_atoms(number, lines, key, cell, lattice_parameter)
            images.append(atoms)

    if index is None:
        return images
    else:
        return images[index]


def make_atoms(index, lines, key, cell, alat=None):
    """Scan through lines to get the atomic positions."""
    atoms = Atoms()
    if key == 'Cartesian axes':
        # initial coordinates are given in terms of 'alat'
        if alat is None:
            warnings.warn("alat expected in make_atoms with cartesian axes")
            alat = 1.0
        for line in lines[index + 3:]:
            entries = line.split()
            if len(entries) == 0:
                break
            symbol = label_to_symbol(entries[1])
            x = float(entries[6])*alat
            y = float(entries[7])*alat
            z = float(entries[8])*alat
            atoms.append(Atom(symbol, (x, y, z)))
        atoms.set_cell(cell)
    elif key == 'ATOMIC_POSITIONS':
        # decide on scale factor based on how the positions
        # are given, choose from:
        # alat, bohr, angstrom, crystal, crystal_sg
        if 'alat' in lines[index]:
            position_scale = alat
            crystal_scale = False
        elif 'bohr' in lines[index]:
            position_scale = units['Bohr']
            crystal_scale = False
        elif 'angstrom' in lines[index]:
            position_scale = 1.0
            crystal_scale = False
        elif 'crystal_sg' in lines[index]:
            raise NotImplementedError('crystal_sg parsing not implemented')
        elif 'crystal' in lines[index]:
            position_scale = 1.0
            crystal_scale = True
        else:
            raise NotImplementedError('{0}'.format(lines[index]))
        for line in lines[index + 1:]:
            entries = line.split()
            if len(entries) == 0 or (entries[0] == 'End'):
                break
            symbol = label_to_symbol(entries[0])
            x = float(entries[1])*position_scale
            y = float(entries[2])*position_scale
            z = float(entries[3])*position_scale
            atoms.append(Atom(symbol, (x, y, z)))
        atoms.set_cell(cell, scale_atoms=crystal_scale)
    # Energy is located after positions.
    energylines = [number for number, line in enumerate(lines) if
                   ('!' in line and 'total energy' in line
                    and number > index)]
    if energylines:
        energyline = min([n for n in energylines if n > index])
        energy = float(lines[energyline].split()[-2]) * units['Ry']
    else:
        energy = None
    # Forces are located after positions.
    forces = np.zeros((len(atoms), 3))
    forcelines = [number for number, line in enumerate(lines) if
                  'Forces acting on atoms (Ry/au):' in line
                  and number > index]
    if forcelines:
        forceline = min([n for n in forcelines if n > index])
        # In QE 5.3 the 'negative rho' has moved above the forceline
        # so need to start 2 lines down.
        if not lines[forceline + 2].strip():
            offset = 4
        else:
            offset = 2
        for line in lines[forceline + offset:]:
            words = line.split()
            if 'force =' in line:
                fx = float(words[-3])
                fy = float(words[-2])
                fz = float(words[-1])
                atom_number = int(words[1]) - 1
                forces[atom_number] = (fx, fy, fz)
            elif len(words) == 0 or 'non-local' in words:
                # 'non-local' line is found with 'high' verbosity
                break
            else:
                continue
        forces *= units['Ry'] / units['Bohr']
    else:
        forces = None
    # Stresses are not always present
    stresslines = [number for number, line in enumerate(lines) if
                   'total   stress  (Ry/bohr**3)' in line and number > index]
    if stresslines:
        stressline = min([n for n in stresslines if n > index])
        xx, xy, xz = lines[stressline + 1].split()[:3]
        yx, yy, yz = lines[stressline + 2].split()[:3]
        zx, zy, zz = lines[stressline + 3].split()[:3]
        stress = np.array([xx, yy, zz, yz, xz, xy], dtype=float)
        stress *= units['Ry'] / (units['Bohr']**3)
    else:
        stress = None
    calc = SinglePointCalculator(atoms, energy=energy, forces=forces,
                                 stress=stress)
    atoms.set_calculator(calc)
    return atoms


def read_espresso_in(fileobj):
    """Parse a Quantum ESPRESSO input files, '.in', '.pwi'.

    ESPRESSO inputs are generally a fortran-namelist format with custom
    blocks of data. The namelist is parsed as a dict and an atoms object
    is constructed from the included information.

    Parameters
    ----------
    fileobj : file | str
        A file-like object that supports line iteration with the contents
        of the input file, or a filename.

    Returns
    -------
    atoms : Atoms
        Structure defined in the input file.

    """
    # TODO: use ase opening mechanisms
    if isinstance(fileobj, basestring):
        fileobj = open(fileobj, 'rU')

    # parse namelist section and extract remaining lines
    data, card_lines = read_fortran_namelist(fileobj)

    # TODO: implemet other ibrav settings
    # get the cell if ibrav=0
    if data['system']['ibrav'] == 0:
        # celldm(1) is in Bohr, A is in angstrom. celldm(1) will be
        # used even if A is also specified.
        if 'celldm(1)' in data['system']:
            alat = data['system']['celldm(1)']*units['Bohr']
        elif 'A' in data['system']:
            alat = data['system']['A']
        else:
            alat = None
        cell = get_cell_parameters(card_lines, alat=alat)
    else:
        raise NotImplementedError('ibrav =/= 0 is not implemented')

    positions_card = get_atomic_positions(
        card_lines, n_atoms=data['system']['nat'], cell=cell, alat=alat)

    symbols = [label_to_symbol(position[0]) for position in positions_card]
    positions = [position[1] for position in positions_card]

    # TODO: put more info into the atoms object
    # e.g magmom, force constraints
    atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)

    return atoms


def get_atomic_positions(lines, n_atoms, cell=None, alat=None):
    """Parse atom positions from ATOMIC_POSITIONS card.

    Parameters
    ----------
    lines : list[str]
        A list of lines containing the ATOMIC_POSITIONS card.
    n_atoms : int
        Expected number of atoms. Only this many lines will be parsed.
    cell : np.array
        Unit cell of the crystal. Only used with crystal coordinates.
    alat : float
        Lattice parameter for atomic coordinated. Only used for alat.

    Returns
    -------
    positions : list[(str, (float, float, float), (float, float, float))]
        A list of the ordered atomic positions in the format:
        label, (x, y, z), (if_x, if_y, if_z)
        Force multipliers are set to None if not present.

    Raises
    ------
    ValueError
        Any problems parsing the data result in ValueError

    """

    positions = None
    # no blanks or comment lines, can the consume n_atoms lines for positions
    trimmed_lines = (line for line in lines
                     if line.strip() and not line[0] == '#')

    for line in trimmed_lines:
        if line.strip().startswith('ATOMIC_POSITIONS'):
            if positions is not None:
                raise ValueError('Multiple ATOMIC_POSITIONS specified')
            # Priority and behaviour tested with QE 5.3
            if 'crystal_sg' in line.lower():
                raise NotImplementedError('CRYSTAL_SG not implemented')
            elif 'crystal' in line.lower():
                cell = cell
            elif 'bohr' in line.lower():
                cell = np.identity(3) * units['Bohr']
            elif 'angstrom' in line.lower():
                cell = np.identity(3)
            #elif 'alat' in line.lower():
            #    cell = np.identity(3) * alat
            else:
                if alat is None:
                    raise ValueError('Set lattice parameter in &SYSTEM for '
                                     'alat coordinates')
                # Always the default, will be DEPRECATED as mandatory
                # in future
                cell = np.identity(3) * alat

            positions = []
            for _atom_idx in range(n_atoms):
                split_line = next(trimmed_lines).split()
                # These can be fractions and other expressions
                position = np.dot((infix_float(split_line[1]),
                                   infix_float(split_line[2]),
                                   infix_float(split_line[3])), cell)
                if len(split_line) > 4:
                    force_mult = (float(split_line[4]),
                                  float(split_line[5]),
                                  float(split_line[6]))
                else:
                    force_mult = None

                positions.append((split_line[0], position, force_mult))

    return positions


def get_cell_parameters(lines, alat=None):
    """Parse unit cell from CELL_PARAMETERS card.

    Parameters
    ----------
    lines : list[str]
        A list with lines containing the CELL_PARAMETERS card.
    alat : float | None
        Unit of lattice vectors in Angstrom. Only used if the card is
        given in units of alat. alat must be None if CELL_PARAMETERS card
        is in Bohr or Angstrom.

    Returns
    -------
    cell : np.array | None
        Cell parameters as a 3x3 array in Angstrom. If no cell is found
        None will be returned instead.

    Raises
    ------
    ValueError
        If CELL_PARAMETERS are given in units of bohr or angstrom
        and alat is not
    """

    cell = None
    # no blanks or comment lines, can take three lines for cell
    trimmed_lines = (line for line in lines
                     if line.strip() and not line[0] == '#')

    for line in trimmed_lines:
        if line.strip().startswith('CELL_PARAMETERS'):
            if cell is not None:
                # multiple definitions
                raise ValueError('CELL_PARAMETERS specified multiple times')
            # Priority and behaviour tested with QE 5.3
            if 'bohr' in line.lower():
                if alat is not None:
                    raise ValueError('Lattice parameters given in '
                                     '&SYSTEM celldm/A and CELL_PARAMETERS '
                                     'bohr')
                cell_units = units['Bohr']
            elif 'angstrom' in line.lower():
                if alat is not None:
                    raise ValueError('Lattice parameters given in '
                                     '&SYSTEM celldm/A and CELL_PARAMETERS '
                                     'angstrom')
                cell_units = 1.0
            elif 'alat' in line.lower():
                if alat is None:
                    raise ValueError('Lattice parameters must be set in '
                                     '&SYSTEM for alat units')
                cell_units = alat
            elif alat is None:
                # may be DEPRECATED in future
                cell_units = units['Bohr']
            else:
                # may be DEPRECATED in future
                cell_units = alat
            # Grab the parameters; blank lines have been removed
            cell = [[ffloat(x) for x in next(trimmed_lines).split()[:3]],
                    [ffloat(x) for x in next(trimmed_lines).split()[:3]],
                    [ffloat(x) for x in next(trimmed_lines).split()[:3]]]
            cell = np.array(cell) * cell_units

    return cell


def str_to_value(string):
    """Attempt to convert string into int, float (including fortran double),
    or bool, in that order, otherwise return the string.
    Valid (case-insensitive) bool values are: '.true.', '.t.', 'true'
    and 't' (or false equivalents).

    Parameters
    ----------
    string : str
        Test to parse for a datatype

    Returns
    -------
    value : any
        Parsed string as the most appropriate datatype of int, float,
        bool or string.

    """

    # Just an integer
    try:
        return int(string)
    except ValueError:
        pass
    # Standard float
    try:
        return float(string)
    except ValueError:
        pass
    # Fortran double
    try:
        return ffloat(string)
    except ValueError:
        pass

    # possible bool, else just the raw string
    if string.lower() in ('.true.', '.t.', 'true', 't'):
        return True
    elif string.lower() in ('.false.', '.f.', 'false', 'f'):
        return False
    else:
        return string.strip("'")


def read_fortran_namelist(fileobj):
    """Takes a fortran-namelist formatted file and returns nested
    dictionaries of sections and key-value data, followed by a list
    of lines of text that do not fit the specifications.

    Behaviour is taken from Quantum ESPRESSO 5.3. Parses fairly
    convoluted files the same way that QE should, but may not get
    all the MANDATORY rules and edge cases for very non-standard files:
        Ignores anything after '!' in a namelist, split pairs on ','
        to include multiple key=values on a line, read values on section
        start and end lines, section terminating character, '/', can appear
        anywhere on a line.
        All of these are ignored if the value is in 'quotes'.

    Parameters
    ----------
    fileobj : file
        An open file-like object.

    Returns
    -------
    data : dict of dict
        Dictionary for each section in the namelist with key = value
        pairs of data.
    card_lines : list of str
        Any lines not used to create the data, assumed to belong to 'cards'
        in the input file.

    """
    # TODO: ignore repeated sections
    # ensure whole file is considered
    fileobj.seek(0)

    # Espresso requires the correct order
    data = OrderedDict()
    card_lines = []
    in_namelist = False
    section = 'none'  # can't be in a section without changing this

    for line in fileobj:
        # leading and trailing whitespace never needed
        line = line.strip()
        if line.startswith('&'):
            # inside a namelist
            section = line.split()[0][1:].lower()  # case insensitive
            data[section] = OrderedDict()
            in_namelist = True
        if not in_namelist and line:
            # TODO, strip comments
            card_lines.append(line)
        if in_namelist:
            # parse k, v from line:
            key = []
            value = None
            in_quotes = False
            for character in line:
                if character == ',' and value is not None and not in_quotes:
                    # finished value:
                    data[section][''.join(key).strip()] = str_to_value(
                        ''.join(value).strip())
                    key = []
                    value = None
                elif character == '=' and value is None and not in_quotes:
                    # start writing value
                    value = []
                elif character == "'":
                    # only found in value anyway
                    in_quotes = not in_quotes
                    value.append("'")
                elif character == '!' and not in_quotes:
                    break
                elif character == '/' and not in_quotes:
                    in_namelist = False
                    break
                elif value is not None:
                    value.append(character)
                else:
                    key.append(character)
            if value is not None:
                data[section][''.join(key).strip()] = str_to_value(
                    ''.join(value).strip())

    return data, card_lines


def ffloat(string):
    """Parse float from fortran compatible float definitions.

    In fortran exponents can be defined with 'd' or 'q' to symbolise
    double or quad precision numbers. Double precision numbers are
    converted to python floats and quad precision values are interpreted
    as numpy longdouble values (platform specific precision).

    Parameters
    ----------
    string : str
        A string containing a number in fortran real format

    Returns
    -------
    value : float | np.longdouble
        Parsed value of the string.

    Raises
    ------
    ValueError
        Unable to parse a float value.

    """

    if 'q' in string.lower():
        return np.longdouble(string.lower().replace('q', 'e'))
    else:
        return float(string.lower().replace('d', 'e'))


def label_to_symbol(label):
    """Convert a valid espresso ATOMIC_SPECIES label to a
    chemical symbol.

    Parameters
    ----------
    label : str
        chemical symbol X (1 or 2 characters, case-insensitive)
        or chemical symbol plus a number or a letter, as in
        "Xn" (e.g. Fe1) or "X_*" or "X-*" (e.g. C1, C_h;
        max total length cannot exceed 3 characters).

    Returns
    -------
    symbol : str
        The best matching species from ase.utils.chemical_symbols

    Raises
    ------
    KeyError
        Couldn't find an appropriate species.

    Notes
    -----
        It's impossible to tell whether e.g. He is helium
        or hydrogen labelled 'e'.
    """

    # possibly a two character species
    # ase Atoms need proper case of chemical symbols.
    if len(label) >= 2:
        test_symbol = label[0].upper() + label[1].lower()
        if test_symbol in chemical_symbols:
            return test_symbol
    # finally try with one character
    test_symbol = label[0].upper()
    if test_symbol in chemical_symbols:
        return test_symbol
    else:
        raise KeyError('Could not parse species from label {0}.'
                       ''.format(label))


def infix_float(text):
    """Parse simple infix maths into a float for compatibility with
    Quantum ESPRESSO ATOMIC_POSITIONS cards. Note: this works with the
    example, and most simple expressions, but the capabilities of
    the two parsers are not identical. Will also parse a normal float
    value properly, but slowly.

    >>> infix_float('1/2*3^(-1/2)')
    0.28867513459481287

    Parameters
    ----------
    text : str
        An arithmetic expression using +, -, *, / and ^, including brackets.

    Returns
    -------
    value : float
        Result of the mathematical expression.

    """

    def middle_brackets(full_text):
        """Extract text from innermost brackets."""
        start, end = 0, len(full_text)
        for (idx, char) in enumerate(full_text):
            if char == '(':
                start = idx
            if char == ')':
                end = idx + 1
                break
        return full_text[start:end]

    def eval_no_bracket_expr(full_text):
        """Calculate value of a mathematical expression, no brackets."""
        exprs = [('+', op.add), ('*', op.mul), ('/', op.div), ('^', op.pow)]
        full_text = full_text.lstrip('(').rstrip(')')
        try:
            return float(full_text)
        except ValueError:
            for symbol, func in exprs:
                if symbol in full_text:
                    left, right = full_text.split(symbol, 1)  # single split
                    return func(eval_no_bracket_expr(left),
                                eval_no_bracket_expr(right))

    while '(' in text:
        middle = middle_brackets(text)
        text = text.replace(middle, '{}'.format(eval_no_bracket_expr(middle)))

    return float(eval_no_bracket_expr(text))

