"""Module to read and write atoms in cif file format.

See http://www.iucr.org/resources/cif/spec/version1.1/cifsyntax for a
description of the file format.  STAR extensions as save frames,
global blocks, nested loops and multi-data values are not supported.
"""

import shlex
import re

import numpy as np

from ase.parallel import paropen
from ase.lattice.spacegroup import crystal
from ase.lattice.spacegroup.spacegroup import spacegroup_from_data



def unread_line(fileobj):
    """Unread the last line read from *fileobj*."""

    # If previous line ends with CRLF, we have to back up one extra 
    # character before entering the loop below
    if fileobj.tell() > 2:
        fileobj.seek(-2, 1)
        if fileobj.read(2) == '\r\n':
            fileobj.seek(-1, 1)

    while True:
        if fileobj.tell() == 0:
            break
        fileobj.seek(-2, 1)
        if fileobj.read(1) in ('\n', '\r'):
            break
        

def convert_value(value):
    """Convert CIF value string to corresponding python type."""
    value = value.strip()
    if re.match('(".*")|(\'.*\')$', value):
        return value[1:-1]
    elif re.match(r'[+-]?\d+$', value):
        return int(value)
    elif re.match(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$', value):
        return float(value)
    elif re.match(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\(\d+\)$', 
                  value):
        return float(value[:value.index('(')])  # strip off uncertainties
    else:
        return value


def parse_singletag(fileobj, line):
    """Parse a CIF tag (entries starting with underscore). Returns
    a key-value pair."""
    kv = line.split(None, 1)
    if len(kv) == 1:  # _key
        key = line
        line = fileobj.readline().strip()
        while line == '':
            line = fileobj.readline().strip()
        if line == ';':  # multiline string 
            lines = [line[1:].lstrip()]
            while True:
                line = fileobj.readline().strip()
                if line == ';':
                    break
                lines.append(line)
            value = '\n'.join(lines).strip()
        else:           # value on next line
            value = line
    else:   # _key, value
        key, value = kv
    return key, convert_value(value)


def parse_loop(fileobj):
    """Parse a CIF loop. Returns a dict with column tag names as keys
    and a lists of the column content as values."""
    header = []
    line = fileobj.readline().strip()
    while line.startswith('_'):
        header.append(line.lower())
        line = fileobj.readline().strip()
    columns = dict([(h, []) for h in header])
    while True:
        lowerline = line.lower()
        if (not line or 
            line.startswith('_') or 
            lowerline.startswith('data_') or 
            lowerline.startswith('loop_')):
            break
        if line.startswith('#'):
            line = fileobj.readline().strip()
            continue
        tokens = shlex.split(line)
        for h, t in zip(header, tokens):
            columns[h].append(convert_value(t))
        line = fileobj.readline().strip()
    if line:
        unread_line(fileobj)
    return columns


def parse_items(fileobj, line):
    """Parse a CIF data items and return a dict with all tags."""
    tags = {}
    while True:
        line = fileobj.readline()
        if not line:
            break
        line = line.strip()
        lowerline = line.lower()
        if not line or line.startswith('#'):
            continue
        elif line.startswith('_'):
            key, value = parse_singletag(fileobj, line)
            tags[key.lower()] = value
        elif lowerline.startswith('loop_'):
            tags.update(parse_loop(fileobj))
        elif lowerline.startswith('data_'):
            unread_line(fileobj)
            break
        else:
            raise ValueError('%s: Unexpected CIF file entry: %s'%(
                    fileobj.name, line))
    return tags


def parse_block(fileobj, line):
    """Parse a CIF data block and return a tuple with the block name
    and a dict with all tags."""
    assert line.lower().startswith('data_')
    blockname = line.split('_', 1)[1].rstrip()
    tags = parse_items(fileobj, line)
    return blockname, tags


def parse_cif(fileobj):
    """Parse a CIF file. Returns a list of blockname and tag
    pairs. All tag names are converted to lower case."""
    if isinstance(fileobj, basestring):
        fileobj = open(fileobj)

    blocks = []
    while True:
        line = fileobj.readline()
        if not line:
            break
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        blocks.append(parse_block(fileobj, line))
    return blocks


def tags2atoms(tags, **kwargs):
    """Returns an Atoms object from a cif tags dictionary."""
    a = tags['_cell_length_a']
    b = tags['_cell_length_b']
    c = tags['_cell_length_c']
    alpha = tags['_cell_angle_alpha']
    beta = tags['_cell_angle_beta']
    gamma = tags['_cell_angle_gamma']

    scaled_positions = np.array([tags['_atom_site_fract_x'], 
                                 tags['_atom_site_fract_y'], 
                                 tags['_atom_site_fract_z']]).T

    symbols = []
    if '_atom_site_type_symbol' in tags:
        labels = tags['_atom_site_type_symbol']
    else:
        labels = tags['_atom_site_label']
    for s in labels:
        # Strip off additional labeling on chemical symbols
        m = re.search(r'([A-Z][a-z]?)', s)  
        symbol = m.group(0)
        symbols.append(symbol)

    # Symmetry specification, see
    # http://www.iucr.org/resources/cif/dictionaries/cif_sym for a
    # complete list of official keys.  In addition we also try to
    # support some commonly used depricated notations
    no = None
    if '_space_group.IT_number' in tags:
        no = tags['_space_group.IT_number']
    elif '_symmetry_int_tables_number' in tags:
        no = tags['_symmetry_int_tables_number']

    symbolHM = None
    if '_space_group.Patterson_name_H-M' in tags:
        symbolHM = tags['_space_group.Patterson_name_H-M']
    elif '_symmetry_space_group_name_H-M' in tags:
        symbolsHM = tags['_symmetry_space_group_name_H-M']

    sitesym = None
    if '_space_group_symop.operation_xyz' in tags:
        sitesym = tags['_space_group_symop.operation_xyz']
    elif '_symmetry_equiv_pos_as_xyz' in tags:
        sitesym = tags['_symmetry_equiv_pos_as_xyz']
        
    spacegroup = 1
    if sitesym is not None:
        spacegroup = spacegroup_from_data(no=no, symbol=symbolHM,
                                          sitesym=sitesym)
    elif no is not None:
        spacegroup = no
    elif symbolHM is not None:
        spacegroup = symbolHM
    else:
        spacegroup = 1

    atoms = crystal(symbols, basis=scaled_positions, 
                    cellpar=[a, b, c, alpha, beta, gamma],
                    spacegroup=spacegroup, **kwargs)
    return atoms
    

def read_cif(fileobj, index=-1, **kwargs):
    """Read Atoms object from CIF file. *index* specifies the data
    block number or name (if string) to return.  Keyword arguments are
    passed on to ase.lattice.spacegroup.crystal()."""
    blocks = parse_cif(fileobj)
    if isinstance(index, str):
        tags = dict(blocks)[index]
        return tags2atoms(tags, **kwargs)
    elif isinstance(index, int):
        name, tags = blocks[index]
        return tags2atoms(tags, **kwargs)
    else:
        return [tags2atoms(tags) for name, tags in blocks[index]]



def write_cif(fileobj, images):
    """Write *images* to CIF file."""
    if isinstance(fileobj, str):
        fileobj = paropen(fileobj, 'w')

    if not isinstance(images, (list, tuple)):
        images = [images]

    for i, atoms in enumerate(images):
        fileobj.write('data_image%d\n' % i)

        from numpy import arccos, pi, dot
        from numpy.linalg import norm

        cell = atoms.cell
        a = norm(cell[0])
        b = norm(cell[1])
        c = norm(cell[2])
        alpha = arccos(dot(cell[1], cell[2])/(b*c))*180./pi
        beta = arccos(dot(cell[0], cell[2])/(a*c))*180./pi
        gamma = arccos(dot(cell[0], cell[1])/(a*b))*180./pi

        fileobj.write('_cell_length_a       %g\n' % a)
        fileobj.write('_cell_length_b       %g\n' % b)
        fileobj.write('_cell_length_c       %g\n' % c)
        fileobj.write('_cell_angle_alpha    %g\n' % alpha)
        fileobj.write('_cell_angle_beta     %g\n' % beta)
        fileobj.write('_cell_angle_gamma    %g\n' % gamma)
        fileobj.write('\n')

        if atoms.pbc.all():
            fileobj.write('_symmetry_space_group_name_H-M    %s\n' % 'P 1')
            fileobj.write('_symmetry_int_tables_number       %d\n' % 1)
            fileobj.write('\n')

            fileobj.write('loop_\n')
            fileobj.write('  _symmetry_equiv_pos_as_xyz\n')
            fileobj.write("  'x, y, z'\n")
            fileobj.write('\n')

        fileobj.write('loop_\n')
        fileobj.write('  _atom_site_label\n')
        fileobj.write('  _atom_site_occupancy\n')
        fileobj.write('  _atom_site_fract_x\n')
        fileobj.write('  _atom_site_fract_y\n')
        fileobj.write('  _atom_site_fract_z\n')
        fileobj.write('  _atom_site_thermal_displace_type\n')
        fileobj.write('  _atom_site_B_iso_or_equiv\n')
        fileobj.write('  _atom_site_type_symbol\n')

        scaled = atoms.get_scaled_positions()
        no = {}
        for i, atom in enumerate(atoms):
            symbol = atom.symbol
            if symbol in no:
                no[symbol] += 1
            else:
                no[symbol] = 1
            fileobj.write(
                '  %-8s %6.4f %7.5f  %7.5f  %7.5f  %4s  %6.3f  %s\n'%(
                    '%s%d' % (symbol, no[symbol]), 
                    1.0, 
                    scaled[i][0], 
                    scaled[i][1], 
                    scaled[i][2],
                    'Biso',
                    1.0,
                    symbol))

    
