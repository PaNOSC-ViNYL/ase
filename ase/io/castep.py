from __future__ import print_function
# -*- coding: utf-8 -*-
"""This module defines I/O routines with CASTEP files.
The key idea is that all function accept or return  atoms objects.
CASTEP specific parameters will be returned through the <atoms>.calc
attribute.
"""
import os
import re
import numpy as np

import ase
from ase.constraints import FixAtoms, FixCartesian
from ase.parallel import paropen

# independent unit management included here:
# When high accuracy is required, this allows to easily pin down
# unit conversion factors from different "unit definition systems"
# (CODATA1986 for ase-3.6.0.2515 vs CODATA2002 for CASTEP 5.01).
#
# ase.units in in ase-3.6.0.2515 is based on CODATA1986
import ase.units
units_ase = {
    'hbar': ase.units._hbar * ase.units.J,
    'Eh': ase.units.Hartree,
    'kB': ase.units.kB,
    'a0': ase.units.Bohr,
    't0': ase.units._hbar * ase.units.J / ase.units.Hartree,
    'c': ase.units._c,
    'me': ase.units._me / ase.units._amu,
    'Pascal': 1.0 / ase.units.Pascal}

# CODATA1986 (included herein for the sake of completeness)
# taken from
#    http://physics.nist.gov/cuu/Archive/1986RMP.pdf
units_CODATA1986 = {
    'hbar': 6.5821220E-16,      # eVs
    'Eh': 27.2113961,           # eV
    'kB': 8.617385E-5,          # eV/K
    'a0': 0.529177249,          # A
    'c': 299792458,             # m/s
    'e': 1.60217733E-19,        # C
    'me': 5.485799110E-4}       # u

# CODATA2002: default in CASTEP 5.01
# (-> check in more recent CASTEP in case of numerical discrepancies?!)
# taken from
#    http://physics.nist.gov/cuu/Document/all_2002.pdf
units_CODATA2002 = {
    'hbar': 6.58211915E-16,     # eVs
    'Eh': 27.2113845,           # eV
    'kB': 8.617343E-5,          # eV/K
    'a0': 0.5291772108,         # A
    'c': 299792458,             # m/s
    'e': 1.60217653E-19,        # C
    'me': 5.4857990945E-4}      # u

# (common) derived entries
for d in (units_CODATA1986, units_CODATA2002):
    d['t0'] = d['hbar'] / d['Eh']     # s
    d['Pascal'] = d['e'] * 1E30       # Pa


__all__ = [
    'read_castep',
    'read_castep_old',
    'read_castep_cell',
    'read_castep_geom',
    'read_md',
    'read_phonon',
    'read_param',
    'read_seed',
    'write_castep_cell',
    'write_param']


def write_castep_cell(filename, atoms, positions_frac=False, castep_cell=None,
                      force_write=False):
    """This CASTEP export function write minimal information to
    a .cell file. If the atoms object is a trajectory, it will
    take the last image.
    """
    if atoms is None:
        print('Atoms object not initialized')
        return False
    if isinstance(atoms, list):
        if len(atoms) > 1:
            atoms = atoms[-1]

    if os.path.isfile(filename) and not force_write:
        print('ase.io.castep.write_param: Set optional argument')
        print('force_write=True to overwrite %s.' % filename)
        return False

    fd = open(filename, 'w')
    fd.write('#######################################################\n')
    fd.write('#CASTEP cell file: %s\n' % filename)
    fd.write('#Created using the Atomic Simulation Environment (ASE)#\n')
    fd.write('#######################################################\n\n')
    fd.write('%BLOCK LATTICE_CART\n')
    cell = np.matrix(atoms.get_cell())
    for line in atoms.get_cell():
        fd.write('    %.10f %.10f %.10f\n' % tuple(line))
    fd.write('%ENDBLOCK LATTICE_CART\n\n\n')

    if positions_frac:
        keyword = 'POSITIONS_FRAC'
        positions = np.array(atoms.get_positions() * cell.I)

    else:
        keyword = 'POSITIONS_ABS'
        positions = atoms.get_positions()

    if (hasattr(atoms, 'calc') and
        hasattr(atoms.calc, 'param') and
        hasattr(atoms.calc.param, 'task')):
        _spin_pol = any([getattr(atoms.calc.param, i).value
                         for i in ['spin_polarized', 'spin_polarised']])
    else:
        _spin_pol = True

    if atoms.get_initial_magnetic_moments().any() and _spin_pol:
        pos_block = [('%s %8.6f %8.6f %8.6f SPIN=%4.2f' %
                     (x, y[0], y[1], y[2], m)) for (x, y, m)
                     in zip(atoms.get_chemical_symbols(),
                            positions,
                            atoms.get_initial_magnetic_moments())]
    else:
        pos_block = [('%s %8.6f %8.6f %8.6f' %
                      (x, y[0], y[1], y[2])) for (x, y)
                     in zip(atoms.get_chemical_symbols(),
                            positions)]

    fd.write('%%BLOCK %s\n' % keyword)
    for line in pos_block:
        fd.write('    %s\n' % line)
    fd.write('%%ENDBLOCK %s\n\n' % keyword)

    # if atoms, has a CASTEP calculator attached, then only
    # write constraints if really necessary
    if (hasattr(atoms, 'calc') and
        hasattr(atoms.calc, 'param') and
        hasattr(atoms.calc.param, 'task')):
        task = atoms.calc.param.task
        if atoms.calc.param.task.value is None:
            suppress_constraints = True
        elif task.value.lower() not in [
            'geometryoptimization',
            # well, CASTEP understands US and UK english...
            'geometryoptimisation',
            'moleculardynamics',
            'transitionstatesearch',
            'phonon']:
            suppress_constraints = True
        else:
            suppress_constraints = False
    else:
        suppress_constraints = True

    constraints = atoms.constraints
    if len(constraints) and not suppress_constraints:
        fd.write('%BLOCK IONIC_CONSTRAINTS \n')
        count = 0
        for constr in constraints:
            if not isinstance(constr, FixAtoms)\
                and not isinstance(constr, FixCartesian)\
                and not suppress_constraints:
                print('Warning: you have constraints in your atoms, that are')
                print('         not supported by CASTEP')
                break
            if isinstance(constr, FixAtoms):
                # sorry, for this complicated block
                # reason is that constraint.index can either
                # hold booleans or integers and in both cases
                # it is an numpy array, so no simple comparison works
                for n, val in enumerate(constr.index):
                    if val.dtype.name.startswith('bool'):
                        if not val:
                            continue
                        symbol = atoms.get_chemical_symbols()[n]
                        nis = atoms.calc._get_number_in_species(n)
                    elif val.dtype.name.startswith('int'):
                        symbol = atoms.get_chemical_symbols()[val]
                        nis = atoms.calc._get_number_in_species(val)
                    else:
                        raise UserWarning('Unrecognized index in' + \
                                           ' constraint %s' % constr)
                    fd.write('%6d %3s %3d   1 0 0 \n' % (count + 1,
                                                         symbol,
                                                         nis))
                    fd.write('%6d %3s %3d   0 1 0 \n' % (count + 2,
                                                         symbol,
                                                         nis))
                    fd.write('%6d %3s %3d   0 0 1 \n' % (count + 3,
                                                         symbol,
                                                         nis))
                    count += 3
            elif isinstance(constr, FixCartesian):
                n = constr.a
                symbol = atoms.get_chemical_symbols()[n]
                nis = atoms.calc._get_number_in_species(n)
                fix_cart = - constr.mask + 1
                if fix_cart[0]:
                    count += 1
                    fd.write('%6d %3s %3d   1 0 0 \n' % (count, symbol, nis))
                if fix_cart[1]:
                    count += 1
                    fd.write('%6d %3s %3d   0 1 0 \n' % (count, symbol, nis))
                if fix_cart[2]:
                    count += 1
                    fd.write('%6d %3s %3d   0 0 1 \n' % (count, symbol, nis))
        fd.write('%ENDBLOCK IONIC_CONSTRAINTS \n')

    if castep_cell is None:
        if hasattr(atoms, 'calc') and hasattr(atoms.calc, 'cell'):
            castep_cell = atoms.calc.cell
        else:
            fd.close()
            return True

    for option in castep_cell._options.values():
        if option.value is not None:
            if option.type == 'Block':
                fd.write('%%BLOCK %s\n' % option.keyword.upper())
                fd.write(option.value)
                fd.write('\n%%ENDBLOCK %s\n' % option.keyword.upper())
            else:
                fd.write('%s : %s\n' % (option.keyword.upper(), option.value))

    fd.close()
    return True


def read_castep_cell(filename, _=None):
    """Read a .cell file and return an atoms object.
    Any value found that does not fit the atoms API
    will be stored in the atoms.calc attribute.

    This routine has been modified to also be able to read *.cell files even if
    there is no CASTEP installation or castep_keywords.py available. We wil
    then make use of a fallback-mode which basically just read atoms positions
    and unit cell information. This can very highly useful for visualization
    using the ASE gui.
    """

    from ase.calculators.castep import Castep

    _fallback = False
    try:
        calc = Castep()
    except Exception as exception:
        print('read_cell: Warning - Was not able to initialize CASTEP calculator.')
        print('           This may be due to a non-existing "castep.keywords.py"')
        print('           file or a non-existing CASTEP installation.')
        print('           Original error message appears below:')
        print('')
        print(' '*11+exception.__str__().replace('\n','\n'+' '*11))
        print('')
        print('           Fallback-mode will be applied to provide at least the')
        print('           geometric information contained in the *.cell file.')
        calc = None
        _fallback = True

    fileobj = open(filename)
    lines = fileobj.readlines()
    fileobj.close()

    def get_tokens(lines, l):
        """Tokenizes one line of a *cell file."""
        comment_chars = '#!'
        while l < len(lines):
            line = lines[l].strip()
            if len(line) == 0:
                l += 1
                continue
            elif any([line.startswith(comment_char)
                      for comment_char in comment_chars]):
                l += 1
                continue
            else:
                for c in comment_chars:
                    if c in line:
                        # icomment = min(line.index(c))
                        # index returns an integer corresponding to the first
                        # appearance, so min() does not work here (and we do
                        # not need it anyways).
                        # see: http://goo.gl/WIEVGs
                        icomment = line.index(c)
                    else:
                        icomment = len(line)
                tokens = line[:icomment].split()
                return tokens, l + 1
        tokens = ''
    # This print statement is definitely not necessary
    #    print("read_cell: Warning - get_tokens has not found any more tokens")
        return tokens, l

    lat = []
    have_lat = False

    pos = []
    spec = []

    # we want to also extract the spins
    magmom = []

    constraints = []
    raw_constraints = {}
    have_pos = False
    pos_frac = False

    l = 0
    while l < len(lines):
        tokens, l = get_tokens(lines, l)
        if not tokens:
            continue
        elif tokens[0].upper() == '%BLOCK':
            if tokens[1].upper() == 'LATTICE_CART' and not have_lat:
                tokens, l = get_tokens(lines, l)
                if len(tokens) == 1:
                    print('read_cell: Warning - ignoring unit specifier in')
                    print('%BLOCK LATTICE_CART (assuming Angstrom instead)')
                    tokens, l = get_tokens(lines, l)
                for _ in range(3):
                    lat_vec = [float(a) for a in tokens[0:3]]
                    lat.append(lat_vec)
                    tokens, l = get_tokens(lines, l)
                if tokens[0].upper() != '%ENDBLOCK':
                    print('read_cell: Warning - ignoring more than three')
                    print('lattice vectors in invalid %BLOCK LATTICE_CART')
                    print('%s ...' % tokens[0].upper())
                have_lat = True

            elif tokens[1].upper() == 'LATTICE_ABC' and not have_lat:
                tokens, l = get_tokens(lines, l)
                if len(tokens) == 1:
                    print('read_cell: Warning - ignoring unit specifier in')
                    print('%BLOCK LATTICE_ABC (assuming Angstrom instead)')
                    tokens, l = get_tokens(lines, l)
                a, b, c = map(float, tokens[0:3])
                tokens, l = get_tokens(lines, l)
                alpha, beta, gamma = [np.radians(float(phi))
                                      for phi in tokens[0:3]]
                tokens, l = get_tokens(lines, l)
                if tokens[0].upper() != '%ENDBLOCK':
                    print('read_cell: Warning - ignoring additional lines in')
                    print('invalid %BLOCK LATTICE_ABC')
                lat_a = [a, 0, 0]
                lat_b = [b * np.cos(gamma), b * np.sin(gamma), 0]
                lat_c1 = c * np.cos(beta)
                lat_c2 = c * ((np.cos(alpha) - np.cos(beta) * np.cos(gamma))
                                 / np.sin(gamma))
                lat_c3 = np.sqrt(c * c - lat_c1 * lat_c1 - lat_c2 * lat_c2)
                lat_c = [lat_c1, lat_c2, lat_c3]
                lat = [lat_a, lat_b, lat_c]
                have_lat = True

            elif tokens[1].upper() == 'POSITIONS_ABS' and not have_pos:
                tokens, l = get_tokens(lines, l)
                if len(tokens) == 1:
                    print('read_cell: Warning - ignoring unit specifier in')
                    print('%BLOCK POSITIONS_ABS(assuming Angstrom instead)')
                    tokens, l = get_tokens(lines, l)
                # fix to be able to read intial spin assigned on the atoms
                while len(tokens) >= 4:
                    spec.append(tokens[0])
                    pos.append([float(p) for p in tokens[1:4]])
                    # read initial spins
                    try:
                        spin = ''.join(tokens[4::]).lower()
                        if not 'spin' in spin:
                            magmom.append(0.)
                        else:
                            magmom.append(float(re.split(r'[:=]+', spin)[-1]))
                    except IndexError:
                        magmom.append(0.)
                    tokens, l = get_tokens(lines, l)
                if tokens[0].upper() != '%ENDBLOCK':
                    print('read_cell: Warning - ignoring invalid lines in')
                    print('%%BLOCK POSITIONS_ABS:\n\t %s' % tokens)
                have_pos = True

            elif tokens[1].upper() == 'POSITIONS_FRAC' and not have_pos:
                pos_frac = True
                tokens, l = get_tokens(lines, l)
                # fix to be able to read intial spin assigned on the atoms
                while len(tokens) >= 4:
                    spec.append(tokens[0])
                    pos.append([float(p) for p in tokens[1:4]])
                    # read initial spins
                    try:
                        spin = ''.join(tokens[4::]).lower()
                        if not 'spin' in spin:
                            magmom.append(0.)
                        else:
                            magmom.append(float(re.split(r'[:=]+', spin)[-1]))
                    except IndexError:
                        magmom.append(0.)
                    tokens, l = get_tokens(lines, l)
                if tokens[0].upper() != '%ENDBLOCK':
                    print('read_cell: Warning - ignoring invalid lines')
                    print('%%BLOCK POSITIONS_FRAC:\n\t %s' % tokens)
                have_pos = True
            elif tokens[1].upper() == 'SPECIES_POT':
                if not _fallback:
                    tokens, l = get_tokens(lines, l)
                    while tokens and not tokens[0].upper() == '%ENDBLOCK':
                        if len(tokens) == 2:
                            calc.cell.species_pot = tuple(tokens)
                        tokens, l = get_tokens(lines, l)
            elif tokens[1].upper() == 'IONIC_CONSTRAINTS':

                while True:
                    if tokens and tokens[0].upper() == '%ENDBLOCK':
                        break
                    tokens, l = get_tokens(lines, l)
                    if not len(tokens) == 6:
                        continue
                    _, species, nic, x, y, z = tokens
                    # convert xyz to floats
                    x = float(x)
                    y = float(y)
                    z = float(z)

                    nic = int(nic)
                    if (species, nic) not in raw_constraints:
                        raw_constraints[(species, nic)] = []
                    raw_constraints[(species, nic)].append(np.array(
                                                           [x, y, z]))

            else:
                print('Warning: the keyword %s is not' % tokens[1].upper())
                print('         interpreted in cell files')
                while not tokens[0].upper() == '%ENDBLOCK':
                    tokens, l = get_tokens(lines, l)
                #raise UserWarning
        else:
            key = tokens[0]
            value = ' '.join(tokens[1:])
            if not _fallback:
                try:
                    calc.__setattr__(key, value)
                except:
                    print('Problem setting calc.cell.%s = %s' % (key, value))
                    raise

    if pos_frac:
        atoms = ase.Atoms(
            calculator=calc,
            cell=lat,
            pbc=True,
            scaled_positions=pos,
            symbols=spec,
            magmoms = magmom
            )
    else:
        atoms = ase.Atoms(
            calculator=calc,
            cell=lat,
            pbc=True,
            positions=pos,
            symbols=spec,
            magmoms = magmom
            )

    fixed_atoms = []
    for (species, nic), value in raw_constraints.items():
        absolute_nr = atoms.calc._get_absolute_number(species, nic)
        if len(value) == 3:
            fixed_atoms.append(absolute_nr)
        elif len(value) == 2:
            constraint = ase.constraints.FixedLine(a=absolute_nr,
                direction=np.cross(value[0], value[1]))
            constraints.append(constraint)
        elif len(value) == 1:
            # catch cases in which constraints are given in a single line in
            # the cell file
            if np.count_nonzero(value[0]) == 3:
                fixed_atoms.append(absolute_nr)
            elif np.count_nonzero(value[0]) == 2:
                # in this case we need a FixedLine instance
                # it is initialized with the atom's index
                constraint = ase.constraints.FixedLine(a=absolute_nr,
                    direction=[not v for v in value[0]])
                constraints.append(constraint)
            else:
                constraint = ase.constraints.FixedPlane(a=absolute_nr,
                    direction=np.array(value[0], dtype=np.float32))
                constraints.append(constraint)
        else:
            print('Error: Found %s statements attached to atoms %s'
                  % (len(value), absolute_nr))
    constraints.append(ase.constraints.FixAtoms(fixed_atoms))
    atoms.set_constraint(constraints)

    if not _fallback:
        # needs to go here again to have the constraints in
        # atoms.calc.atoms.constraints as well
        atoms.calc.atoms = atoms
        atoms.calc.push_oldstate()
    return atoms

def read_castep(filename, _ = None):
    """
    This routine replaces the former read_castep() routine. Basically it does
    the same job, but it uses the read() functionality from the Castep
    calculator class. This allows a much more complete parsing and we do not
    have to take care of syncing the respective routine with each other.

    Note: This routine returns a single atoms_object only, whereas the former
    routine, in principle, returned a list of atoms objects. Yet, if you want
    to parse an MD run, use the novel function `read_md()`

    There is no use of the "index" argument as of now, it is just inserted for
    convenience to comply with the generic "read()" in ase.io

    Please note that this routine will return an atom ordering as found
    within the castep file. This means that the species will be ordered by
    ascending atomic numbers. The atoms witin a species are ordered as given
    in the original cell file.
    """
    from ase.calculators.castep import Castep

    calc = Castep()
    calc.read(castep_file=filename)

    # now we trick the calculator instance such that we can savely extract
    # energies and forces from this atom. Basically what we do is to trick the
    # internal routine calculation_required() to always return False such that
    # we do not need to re-run a CASTEP calculation.
    #
    # Probably we can solve this with a flag to the read() routine at some
    # point, but for the moment I do not want to change too much in there.
    calc._old_atoms = calc.atoms
    calc._old_param = calc.param
    calc._old_cell = calc.cell

    return calc.atoms


# this actually does not belong here
# think how one could join this with
# the ase.calculators.castep.Castep.read()
# in the future!
# --> has been done (see above)

def read_castep_old(filename, _=-1):
    """Reads a .castep file and returns an atoms  object.
    The calculator information will be stored in the calc attribute.
    If more than one SCF step is found, a list of all steps
    will be stored in the traj attribute.

    Note that the index argument has no effect as of now.

    Please note that this routine will return an atom ordering as found
    within the castep file. This means that the species will be ordered by
    ascending atomic numbers. The atoms witin a species are ordered as given
    in the original cell file.
    """
    from ase.calculators.singlepoint import SinglePointCalculator

    fileobj = open(filename)
    lines = fileobj.readlines()
    fileobj.close()
    traj = []
    energy_total = None
    energy_0K = None
    for i, line in enumerate(lines):
        if 'NB est. 0K energy' in line:
            energy_0K = float(line.split()[6])
        # support also for dispersion correction
        elif 'NB dispersion corrected est. 0K energy*' in line:
            energy_0K = float(line.split()[-2])
        elif 'Final energy, E' in line:
            energy_total = float(line.split()[4])
        elif 'Dispersion corrected final energy' in line:
            pass
            # dispcorr_energy_total = float(line.split()[-2])
            # sedc_apply = True
        elif 'Dispersion corrected final free energy' in line:
            pass  # dispcorr_energy_free = float(line.split()[-2])
        elif 'dispersion corrected est. 0K energy' in line:
            pass  # dispcorr_energy_0K = float(line.split()[-2])
        elif 'Unit Cell' in line:
            cell = [x.split()[0:3] for x in lines[i + 3:i + 6]]
            cell = np.array([[float(col) for col in row] for row in cell])
        elif 'Cell Contents' in line:
            geom_starts = i
            start_found = False
            for j, jline in enumerate(lines[geom_starts:]):
                if jline.find('xxxxx') > 0 and start_found:
                    geom_stop = j + geom_starts
                    break
                if jline.find('xxxx') > 0 and not start_found:
                    geom_start = j + geom_starts + 4
                    start_found = True
            species = [line.split()[1] for line in lines[geom_start:geom_stop]]
            geom = np.dot(np.array([[float(col) for col in line.split()[3:6]]
                                    for line in lines[geom_start:geom_stop]]),
                          cell)
        elif 'Writing model to' in line:
            atoms = ase.Atoms(
                cell=cell,
                pbc=True,
                positions=geom,
                symbols=''.join(species))
            # take 0K energy where available, else total energy
            if energy_0K:
                energy = energy_0K
            else:
                energy = energy_total
            # generate a minimal single-point calculator
            sp_calc = SinglePointCalculator(atoms=atoms,
                                            energy=energy,
                                            forces=None,
                                            magmoms=None,
                                            stress=None)
            atoms.set_calculator(sp_calc)
            traj.append(atoms)

    return traj


def read_param(filename, calc=None):
    """Reads a param file. If an Castep object is passed as the
    second argument, the parameter setings are merged into
    the existing object and returned. Otherwise a new Castep()
    calculator instance gets created and returned.

    Parameters:
        filename: the .param file. Only opens reading
        calc: [Optional] calculator object to hang parameters onto
    """
    if calc is None:
        from ase.calculators.castep import Castep
        calc = Castep(check_castep_version=False)
    calc.merge_param(filename)
    return calc


def write_param(filename, param, check_checkfile=False,
                                 force_write=False,
                                 interface_options=None):
    """Writes a CastepParam object to a CASTEP .param file

    Parameters:
        filename: the location of the file to write to. If it
        exists it will be overwritten without warning. If it
        doesn't it will be created.
        param: a CastepParam instance
        check_checkfile : if set to True, write_param will
        only write continuation or reuse statement
        if a restart file exists in the same directory
    """
    if os.path.isfile(filename) and not force_write:
        print('ase.io.castep.write_param: Set optional argument')
        print('force_write=True to overwrite %s.' % filename)
        return False

    out = paropen(filename, 'w')
    out.write('#######################################################\n')
    out.write('#CASTEP param file: %s\n' % filename)
    out.write('#Created using the Atomic Simulation Environment (ASE)#\n')
    if interface_options is not None:
        out.write('# Internal settings of the calculator\n')
        out.write('# This can be switched off by settings\n')
        out.write('# calc._export_settings = False\n')
        out.write('# If stated, this will be automatically processed\n')
        out.write('# by ase.io.castep.read_seed()\n')
        for option, value in sorted(interface_options.items()):
            out.write('# ASE_INTERFACE %s : %s\n' % (option, value))
    out.write('#######################################################\n\n')
    for keyword, opt in sorted(param._options.items()):
        if opt.type == 'Defined':
            if opt.value is not None:
                out.write('%s\n' % (option))
        elif opt.value is not None:
            if keyword in ['continuation', 'reuse'] and check_checkfile:
                if opt.value == 'default':
                    if not os.path.exists('%s.%s'\
                        % (os.path.splitext(filename)[0], 'check')):
                        continue
                elif not (os.path.exists(opt.value)
                          # CASTEP also understands relative path names, hence
                          # also check relative to the param file directory
                          or os.path.exists(os.path.join(
                                                os.path.dirname(filename),
                                                opt.value))
                          ):
                    continue
            out.write('%s : %s\n'
                % (keyword, opt.value))
    out.close()


def read_castep_geom(filename, _=-1, units=units_CODATA2002):
    """Reads a .geom file produced by the CASTEP GeometryOptimization task and
    returns an atoms  object.
    The information about total free energy and forces of each atom for every
    relaxation step will be stored for further analysis especially in a
    single-point calculator.
    Note that everything in the .geom file is in atomic units, which has
    been conversed to commonly used unit angstrom(length) and eV (energy).

    Note that the index argument has no effect as of now.

    Contribution by Wei-Bing Zhang. Thanks!
    """
    from ase.calculators.singlepoint import SinglePointCalculator

    # extended support for compressed formats
    fname = filename.lower()
    if fname.endswith('.gz'):
        import gzip
        fileobj = gzip.open(filename)
    elif fname.endswith('.bz2'):
        import bz2
        fileobj = bz2.BZ2File(filename)
    else:
        fileobj = open(filename)
    txt = fileobj.readlines()
    fileobj.close()

    traj = []

    Hartree = units['Eh']
    Bohr = units['a0']

    # Yeah, we know that...
    #print('N.B.: Energy in .geom file is not 0K extrapolated.')
    for i, line in enumerate(txt):
        if line.find('<-- E') > 0:
            start_found = True
            energy = float(line.split()[0]) * Hartree
            cell = [x.split()[0:3] for x in txt[i + 1:i + 4]]
            cell = np.array([[float(col) * Bohr for col in row] for row in
                cell])
        if line.find('<-- R') > 0 and start_found:
            start_found = False
            geom_start = i
            for i, line in enumerate(txt[geom_start:]):
                if line.find('<-- F') > 0:
                    geom_stop = i + geom_start
                    break
            species = [line.split()[0] for line in
                txt[geom_start:geom_stop]]
            geom = np.array([[float(col) * Bohr for col in
                line.split()[2:5]] for line in txt[geom_start:geom_stop]])
            forces = np.array([[float(col) * Hartree / Bohr for col in
                line.split()[2:5]] for line in
                    txt[geom_stop:geom_stop + (geom_stop - geom_start)]])
            image = ase.Atoms(species, geom, cell=cell, pbc=True)
            image.set_calculator(SinglePointCalculator(energy, forces, None,
                None, image))
            traj.append(image)

    return traj


def read_md(filename, _=-1, return_scalars=False, units=units_CODATA2002):
    """Reads a .md file written by a CASTEP MolecularDynamics task
    and returns the trajectory stored therein as a list of atoms object.

    Note that the index argument has no effect as of now."""

    from ase.calculators.singlepoint import SinglePointCalculator

    factors = {
               't': units['t0'] * 1E15,     # fs
               'E': units['Eh'],            # eV
               'T': units['Eh'] / units['kB'],
               'P': units['Eh'] / units['a0']**3 * units['Pascal'],
               'h': units['a0'],
               'hv':units['a0'] / units['t0'],
               'S': units['Eh'] / units['a0']**3,
               'R': units['a0'],
               'V': np.sqrt(units['Eh'] / units['me']),
               'F': units['Eh'] / units['a0']
               }

    fname = filename.lower()
    if fname.endswith('.gz'):
        import gzip
        f = gzip.open(filename)
    elif fname.endswith('.bz2'):
        import bz2
        f = bz2.BZ2File(filename)
    else:
        f = open(filename)
    lines = f.readlines()
    f.close()

    l = 0
    while 'END header' not in lines[l]:
        l += 1
    l_end_header = l
    lines = lines[l_end_header + 1:]
    times = []
    energies = []
    temperatures = []
    pressures = []
    traj = []

    # Initialization
    time = None
    Epot = None
    Ekin = None
    EH = None
    temperature = None
    pressure = None
    symbols = None
    positions = None
    cell = None
    velocities = None
    symbols = []
    positions = []
    velocities = []
    forces = []
    cell = np.eye(3)
    cell_velocities = []
    stress = []

    for (l, line) in enumerate(lines):
        fields = line.split()
        if len(fields) == 0:
            if l != 0:
                times.append(time)
                energies.append([Epot, EH, Ekin])
                temperatures.append(temperature)
                pressures.append(pressure)
                atoms = ase.Atoms(symbols=symbols,
                                  positions=positions,
                                  cell=cell)
                atoms.set_velocities(velocities)
                if len(stress) == 0:
                    atoms.set_calculator(SinglePointCalculator(Epot,
                                                forces, None, None, atoms))
                else:
                    atoms.set_calculator(SinglePointCalculator(Epot,
                                                forces, stress, None, atoms))
                traj.append(atoms)
            symbols = []
            positions = []
            velocities = []
            forces = []
            cell = []
            cell_velocities = []
            stress = []
            continue
        if len(fields) == 1:
            time = factors['t'] * float(fields[0])
            continue

        if fields[-1] == 'E':
            E = [float(x) for x in fields[0:3]]
            Epot, EH, Ekin = [factors['E'] * Ei for Ei in E]
            continue

        if fields[-1] == 'T':
            temperature = factors['T'] * float(fields[0])
            continue

        # only printed in case of variable cell calculation or calculate_stress
        # explicitly requested
        if fields[-1] == 'P':
            pressure = factors['P'] * float(fields[0])
            continue
        if fields[-1] == 'h':
            h = [float(x) for x in fields[0:3]]
            cell.append([factors['h'] * hi for hi in h])
            continue

        # only printed in case of variable cell calculation
        if fields[-1] == 'hv':
            hv = [float(x) for x in fields[0:3]]
            cell_velocities.append([factors['hv'] * hvi for hvi in hv])
            continue

        # only printed in case of variable cell calculation
        if fields[-1] == 'S':
            S = [float(x) for x in fields[0:3]]
            stress.append([factors['S'] * Si for Si in S])
            continue
        if fields[-1] == 'R':
            symbols.append(fields[0])
            R = [float(x) for x in fields[2:5]]
            positions.append([factors['R'] * Ri for Ri in R])
            continue
        if fields[-1] == 'V':
            V = [float(x) for x in fields[2:5]]
            velocities.append([factors['V'] * Vi for Vi in V])
            continue
        if fields[-1] == 'F':
            F = [float(x) for x in fields[2:5]]
            forces.append([factors['F'] * Fi for Fi in F])
            continue

    if return_scalars:
        data = [times, energies, temperatures, pressures]
        return data, traj
    else:
        return traj


def read_phonon(filename, _=-1, read_vib_data=False,
                gamma_only=True, frequency_factor=None,
                units=units_CODATA2002):
    """
    Reads a .phonon file written by a CASTEP Phonon task and returns an atoms
    object, as well as the calculated vibrational data if requested.

    Note that the index argument has no effect as of now."""

    fname = filename.lower()
    if fname.endswith('.gz'):
        import gzip
        f = gzip.open(filename)
    elif fname.endswith('.bz2'):
        import bz2
        f = bz2.BZ2File(filename)
    else:
        f = open(filename)
    lines = f.readlines()
    f.close()

    atoms = None
    cell = []
    N = Nb = Nq = 0
    scaled_positions = []
    symbols = []
    masses = []

    # header
    l = 0
    while l < len(lines):

        line = lines[l]

        if 'Number of ions' in line:
            N = int(line.split()[3])
        elif 'Number of branches' in line:
            Nb = int(line.split()[3])
        elif 'Number of wavevectors'in line:
            Nq = int(line.split()[3])
        elif 'Unit cell vectors (A)' in line:
            for ll in range(3):
                l += 1
                fields = lines[l].split()
                cell.append([float(x) for x in fields[0:3]])
        elif 'Fractional Co-ordinates' in line:
            for ll in range(N):
                l += 1
                fields = lines[l].split()
                scaled_positions.append([float(x) for x in fields[1:4]])
                symbols.append(fields[4])
                masses.append(float(fields[5]))
        elif 'END header' in line:
            l += 1
            atoms = ase.Atoms(symbols=symbols,
                              scaled_positions=scaled_positions,
                              cell=cell)
            break

        l += 1

    # Eigenmodes and -vectors
    if frequency_factor is None:
        Kayser_to_eV = 1E2 * 2 * np.pi * units['hbar'] * units['c']
    # N.B. "fixed default" unit for frequencies in .phonon files is "cm-1"
    # (i.e. the latter is unaffected by the internal unit conversion system of
    # CASTEP!) set conversion factor to convert therefrom to eV by default for
    # now
    frequency_factor = Kayser_to_eV
    qpoints = []
    weights = []
    frequencies = []
    displacements = []
    for nq in range(Nq):
        fields = lines[l].split()
        qpoints.append([float(x) for x in fields[2:5]])
        weights.append(float(fields[5]))
    freqs = []
    for ll in range(Nb):
        l += 1
        fields = lines[l].split()
        freqs.append(frequency_factor * float(fields[1]))
    frequencies.append(np.array(freqs))

    # skip the two Phonon Eigenvectors header lines
    l += 2

    # generate a list of displacements with a structure that is identical to
    # what is stored internally in the Vibrations class (see in
    # ase.vibrations.Vibrations.modes):
    #      np.array(displacements).shape == (Nb,3*N)

    disps = []
    for ll in range(Nb):
        disp_coords = []
        for lll in range(N):
            l += 1
            fields = lines[l].split()
            disp_x = float(fields[2]) + float(fields[3]) * 1.0j
            disp_y = float(fields[4]) + float(fields[5]) * 1.0j
            disp_z = float(fields[6]) + float(fields[7]) * 1.0j
            disp_coords.extend([disp_x, disp_y, disp_z])
        disps.append(np.array(disp_coords))
    displacements.append(np.array(disps))

    if read_vib_data:
        if gamma_only:
            vibdata = [frequencies[0], displacements[0]]
        else:
            vibdata = [qpoints, weights, frequencies, displacements]
        return vibdata, atoms
    else:
        return atoms


def read_seed(seed, new_seed=None, ignore_internal_keys=False):
    """A wrapper around the CASTEP Calculator in conjunction with
    read_cell and read_param. Basically this can be used to reuse
    a previous calculation which results in a triple of
    cell/param/castep file. The label of the calculation if pre-
    fixed with `copy_of_` and everything else will be recycled as
    much as possible from the addressed calculation.

    Please note that this routine will return an atoms ordering as specified
    in the cell file! It will thus undo the potential reordering internally
    done by castep.
    """

    directory = os.path.abspath(os.path.dirname(seed))
    seed = os.path.basename(seed)

    paramfile = os.path.join(directory, '%s.param' % seed)
    cellfile = os.path.join(directory, '%s.cell' % seed)
    castepfile = os.path.join(directory, '%s.castep' % seed)

    atoms = read_castep_cell(cellfile)
    atoms.calc._directory = directory
    atoms.calc._rename_existing_dir = False
    atoms.calc._castep_pp_path = directory
    atoms.calc.merge_param(paramfile,
                           ignore_internal_keys=ignore_internal_keys)
    if new_seed is None:
        atoms.calc._label = 'copy_of_%s' % seed
    else:
        atoms.calc._label = str(new_seed)
    if os.path.isfile(castepfile):
        # _set_atoms needs to be True here
        # but we set it right back to False
        #atoms.calc._set_atoms = False
        # BUGFIX: I do not see a reason to do that!
        atoms.calc.read(castepfile)
        #atoms.calc._set_atoms = False

        # sync the top-level object with the
        # one attached to the calculator
        atoms = atoms.calc.atoms
    else:
        # There are cases where we only want to restore a calculator/atoms
        # setting without a castep file...
        pass
        # No print statement required in these cases
        print('Corresponding *.castep file not found.')
        print('Atoms object will be restored from *.cell and *.param only.')
    atoms.calc.push_oldstate()

    return atoms
