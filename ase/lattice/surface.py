"""Helper functions for creating the most common surfaces and related tasks.

The helper functions can create the most common low-index surfaces,
add vacuum layers and add adsorbates.

"""

from __future__ import division
from math import sqrt, hypot, acos
from operator import itemgetter

import numpy as np

from ase.atom import Atom
from ase.atoms import Atoms
from ase.data import reference_states, atomic_numbers
from ase.lattice.cubic import FaceCenteredCubic
from ase.lattice.general_surface import surface

__all__ = ['surface', 'add_adsorbate', 'add_vacuum',
           'bcc100', 'bcc110', 'bcc111',
           'diamond100', 'diamond111',
           'fcc100', 'fcc110', 'fcc111', 'fcc211',
           'hcp0001', 'hcp10m10', 'mx2']


def fcc100(symbol, size, a=None, vacuum=None):
    """FCC(100) surface.
 
    Supported special adsorption sites: 'ontop', 'bridge', 'hollow'."""
    return _surface(symbol, 'fcc', '100', size, a, None, vacuum)


def fcc110(symbol, size, a=None, vacuum=None):
    """FCC(110) surface.
 
    Supported special adsorption sites: 'ontop', 'longbridge',
    'shortbridge','hollow'."""
    return _surface(symbol, 'fcc', '110', size, a, None, vacuum)


def bcc100(symbol, size, a=None, vacuum=None):
    """BCC(100) surface.
 
    Supported special adsorption sites: 'ontop', 'bridge', 'hollow'."""
    return _surface(symbol, 'bcc', '100', size, a, None, vacuum)


def bcc110(symbol, size, a=None, vacuum=None, orthogonal=False):
    """BCC(110) surface.
 
    Supported special adsorption sites: 'ontop', 'longbridge',
    'shortbridge', 'hollow'.
 
    Use *orthogonal=True* to get an orthogonal unit cell - works only
    for size=(i,j,k) with j even."""
    return _surface(symbol, 'bcc', '110', size, a, None, vacuum, orthogonal)


def bcc111(symbol, size, a=None, vacuum=None, orthogonal=False):
    """BCC(111) surface.
 
    Supported special adsorption sites: 'ontop'.
 
    Use *orthogonal=True* to get an orthogonal unit cell - works only
    for size=(i,j,k) with j even."""
    return _surface(symbol, 'bcc', '111', size, a, None, vacuum, orthogonal)


def fcc111(symbol, size, a=None, vacuum=None, orthogonal=False):
    """FCC(111) surface.
 
    Supported special adsorption sites: 'ontop', 'bridge', 'fcc' and 'hcp'.
 
    Use *orthogonal=True* to get an orthogonal unit cell - works only
    for size=(i,j,k) with j even."""
    return _surface(symbol, 'fcc', '111', size, a, None, vacuum, orthogonal)


def hcp0001(symbol, size, a=None, c=None, vacuum=None, orthogonal=False):
    """HCP(0001) surface.
 
    Supported special adsorption sites: 'ontop', 'bridge', 'fcc' and 'hcp'.
 
    Use *orthogonal=True* to get an orthogonal unit cell - works only
    for size=(i,j,k) with j even."""
    return _surface(symbol, 'hcp', '0001', size, a, c, vacuum, orthogonal)

    
def hcp10m10(symbol, size, a=None, c=None, vacuum=None):
    """HCP(10m10) surface.
    
    Supported special adsorption sites: 'ontop'.
    
    Works only for size=(i,j,k) with j even."""
    return _surface(symbol, 'hcp', '10m10', size, a, c, vacuum)


def diamond100(symbol, size, a=None, vacuum=None):
    """DIAMOND(100) surface.

    Supported special adsorption sites: 'ontop'."""
    return _surface(symbol, 'diamond', '100', size, a, None, vacuum)


def diamond111(symbol, size, a=None, vacuum=None, orthogonal=False):
    """DIAMOND(111) surface.
 
    Supported special adsorption sites: 'ontop'."""

    if orthogonal:
        raise NotImplementedError("Can't do orthogonal cell yet!")
    return _surface(symbol, 'diamond', '111', size, a, None, vacuum,
                    orthogonal)

    
def add_adsorbate(slab, adsorbate, height, position=(0, 0), offset=None,
                  mol_index=0):
    """Add an adsorbate to a surface.

    This function adds an adsorbate to a slab.  If the slab is
    produced by one of the utility functions in ase.lattice.surface, it
    is possible to specify the position of the adsorbate by a keyword
    (the supported keywords depend on which function was used to
    create the slab).

    If the adsorbate is a molecule, the atom indexed by the mol_index
    optional argument is positioned on top of the adsorption position
    on the surface, and it is the responsibility of the user to orient
    the adsorbate in a sensible way.

    This function can be called multiple times to add more than one
    adsorbate.

    Parameters:

    slab: The surface onto which the adsorbate should be added.

    adsorbate:  The adsorbate. Must be one of the following three types:
        A string containing the chemical symbol for a single atom.
        An atom object.
        An atoms object (for a molecular adsorbate).

    height: Height above the surface.

    position: The x-y position of the adsorbate, either as a tuple of
        two numbers or as a keyword (if the surface is produced by one
        of the functions in ase.lattice.surfaces).

    offset (default: None): Offsets the adsorbate by a number of unit
        cells. Mostly useful when adding more than one adsorbate.

    mol_index (default: 0): If the adsorbate is a molecule, index of
        the atom to be positioned above the location specified by the
        position argument.

    Note *position* is given in absolute xy coordinates (or as
    a keyword), whereas offset is specified in unit cells.  This
    can be used to give the positions in units of the unit cell by
    using *offset* instead.
    
    """
    info = slab.adsorbate_info
    if 'cell' not in info:
        info['cell'] = slab.get_cell()[:2, :2]

    pos = np.array([0.0, 0.0])  # (x, y) part
    spos = np.array([0.0, 0.0])  # part relative to unit cell
    if offset is not None:
        spos += np.asarray(offset, float)

    if isinstance(position, str):
        # A site-name:
        if 'sites' not in info:
            raise TypeError('If the atoms are not made by an ' +
                            'ase.lattice.surface function, ' +
                            'position cannot be a name.')
        if position not in info['sites']:
            raise TypeError('Adsorption site %s not supported.' % position)
        spos += info['sites'][position]
    else:
        pos += position

    pos += np.dot(spos, info['cell'])

    # Convert the adsorbate to an Atoms object
    if isinstance(adsorbate, Atoms):
        ads = adsorbate
    elif isinstance(adsorbate, Atom):
        ads = Atoms([adsorbate])
    else:
        # Assume it is a string representing a single Atom
        ads = Atoms([Atom(adsorbate)])

    # Get the z-coordinate:
    try:
        a = info['top layer atom index']
    except KeyError:
        a = slab.positions[:, 2].argmax()
        info['top layer atom index'] = a
    z = slab.positions[a, 2] + height

    # Move adsorbate into position
    ads.translate([pos[0], pos[1], z] - ads.positions[mol_index])

    # Attach the adsorbate
    slab.extend(ads)

    
def add_vacuum(atoms, vacuum):
    """Add vacuum layer to the atoms.

    Parameters:
 
    atoms: An Atoms object most likely created by one of the
    ase.lattice modules.
 
    vacuum: The thickness of the vacuum layer (in Angstrom).
    """
    uc = atoms.get_cell()
    normal = np.cross(uc[0], uc[1])
    costheta = np.dot(normal, uc[2]) / np.sqrt(np.dot(normal, normal) *
                                               np.dot(uc[2], uc[2]))
    length = np.sqrt(np.dot(uc[2], uc[2]))
    newlength = length + vacuum / costheta
    uc[2] *= newlength / length
    atoms.set_cell(uc)
    
    
def _surface(symbol, structure, face, size, a, c, vacuum, orthogonal=True):
    """Function to build often used surfaces.

    Don't call this function directly - use fcc100, fcc110, bcc111, ..."""
    
    Z = atomic_numbers[symbol]

    if a is None:
        sym = reference_states[Z]['symmetry']
        if sym != structure:
            raise ValueError("Can't guess lattice constant for %s-%s!" %
                             (structure, symbol))
        a = reference_states[Z]['a']

    if structure == 'hcp' and c is None:
        if reference_states[Z]['symmetry'] == 'hcp':
            c = reference_states[Z]['c/a'] * a
        else:
            c = sqrt(8 / 3.0) * a

    positions = np.empty((size[2], size[1], size[0], 3))
    positions[..., 0] = np.arange(size[0]).reshape((1, 1, -1))
    positions[..., 1] = np.arange(size[1]).reshape((1, -1, 1))
    positions[..., 2] = np.arange(size[2]).reshape((-1, 1, 1))

    numbers = np.ones(size[0] * size[1] * size[2], int) * Z

    tags = np.empty((size[2], size[1], size[0]), int)
    tags[:] = np.arange(size[2], 0, -1).reshape((-1, 1, 1))

    slab = Atoms(numbers,
                 tags=tags.ravel(),
                 pbc=(True, True, False),
                 cell=size)

    surface_cell = None
    sites = {'ontop': (0, 0)}
    surf = structure + face
    if surf == 'fcc100':
        cell = (sqrt(0.5), sqrt(0.5), 0.5)
        positions[-2::-2, ..., :2] += 0.5
        sites.update({'hollow': (0.5, 0.5), 'bridge': (0.5, 0)})
    elif surf == 'diamond100':
        cell = (sqrt(0.5), sqrt(0.5), 0.5 / 2)
        positions[-4::-4, ..., :2] += (0.5, 0.5)
        positions[-3::-4, ..., :2] += (0.0, 0.5)
        positions[-2::-4, ..., :2] += (0.0, 0.0)
        positions[-1::-4, ..., :2] += (0.5, 0.0)
    elif surf == 'fcc110':
        cell = (1.0, sqrt(0.5), sqrt(0.125))
        positions[-2::-2, ..., :2] += 0.5
        sites.update({'hollow': (0.5, 0.5), 'longbridge': (0.5, 0),
                      'shortbridge': (0, 0.5)})
    elif surf == 'bcc100':
        cell = (1.0, 1.0, 0.5)
        positions[-2::-2, ..., :2] += 0.5
        sites.update({'hollow': (0.5, 0.5), 'bridge': (0.5, 0)})
    else:
        if orthogonal and size[1] % 2 == 1:
            raise ValueError(("Can't make orthorhombic cell with size=%r.  " %
                              (tuple(size),)) +
                             'Second number in size must be even.')
        if surf == 'fcc111':
            cell = (sqrt(0.5), sqrt(0.375), 1 / sqrt(3))
            if orthogonal:
                positions[-1::-3, 1::2, :, 0] += 0.5
                positions[-2::-3, 1::2, :, 0] += 0.5
                positions[-3::-3, 1::2, :, 0] -= 0.5
                positions[-2::-3, ..., :2] += (0.0, 2.0 / 3)
                positions[-3::-3, ..., :2] += (0.5, 1.0 / 3)
            else:
                positions[-2::-3, ..., :2] += (-1.0 / 3, 2.0 / 3)
                positions[-3::-3, ..., :2] += (1.0 / 3, 1.0 / 3)
            sites.update({'bridge': (0.5, 0), 'fcc': (1.0 / 3, 1.0 / 3),
                          'hcp': (2.0 / 3, 2.0 / 3)})
        elif surf == 'diamond111':
            cell = (sqrt(0.5), sqrt(0.375), 1 / sqrt(3) / 2)
            assert not orthogonal
            positions[-1::-6, ..., :3] += (0.0, 0.0, 0.5)
            positions[-2::-6, ..., :2] += (0.0, 0.0)
            positions[-3::-6, ..., :3] += (-1.0 / 3, 2.0 / 3, 0.5)
            positions[-4::-6, ..., :2] += (-1.0 / 3, 2.0 / 3)
            positions[-5::-6, ..., :3] += (1.0 / 3, 1.0 / 3, 0.5)
            positions[-6::-6, ..., :2] += (1.0 / 3, 1.0 / 3)
        elif surf == 'hcp0001':
            cell = (1.0, sqrt(0.75), 0.5 * c / a)
            if orthogonal:
                positions[:, 1::2, :, 0] += 0.5
                positions[-2::-2, ..., :2] += (0.0, 2.0 / 3)
            else:
                positions[-2::-2, ..., :2] += (-1.0 / 3, 2.0 / 3)
            sites.update({'bridge': (0.5, 0), 'fcc': (1.0 / 3, 1.0 / 3),
                          'hcp': (2.0 / 3, 2.0 / 3)})
        elif surf == 'hcp10m10':
            cell = (1.0, 0.5 * c / a, sqrt(0.75))
            assert orthogonal
            positions[-2::-2, ..., 0] += 0.5
            positions[:, ::2, :, 2] += 2.0 / 3
        elif surf == 'bcc110':
            cell = (1.0, sqrt(0.5), sqrt(0.5))
            if orthogonal:
                positions[:, 1::2, :, 0] += 0.5
                positions[-2::-2, ..., :2] += (0.0, 1.0)
            else:
                positions[-2::-2, ..., :2] += (-0.5, 1.0)
            sites.update({'shortbridge': (0, 0.5),
                          'longbridge': (0.5, 0),
                          'hollow': (0.375, 0.25)})
        elif surf == 'bcc111':
            cell = (sqrt(2), sqrt(1.5), sqrt(3) / 6)
            if orthogonal:
                positions[-1::-3, 1::2, :, 0] += 0.5
                positions[-2::-3, 1::2, :, 0] += 0.5
                positions[-3::-3, 1::2, :, 0] -= 0.5
                positions[-2::-3, ..., :2] += (0.0, 2.0 / 3)
                positions[-3::-3, ..., :2] += (0.5, 1.0 / 3)
            else:
                positions[-2::-3, ..., :2] += (-1.0 / 3, 2.0 / 3)
                positions[-3::-3, ..., :2] += (1.0 / 3, 1.0 / 3)
            sites.update({'hollow': (1.0 / 3, 1.0 / 3)})
        else:
            2 / 0
            
        surface_cell = a * np.array([(cell[0], 0),
                                     (cell[0] / 2, cell[1])])
        if not orthogonal:
            cell = np.array([(cell[0], 0, 0),
                             (cell[0] / 2, cell[1], 0),
                             (0, 0, cell[2])])

    if surface_cell is None:
        surface_cell = a * np.diag(cell[:2])

    if isinstance(cell, tuple):
        cell = np.diag(cell)
        
    slab.set_positions(positions.reshape((-1, 3)))

    slab.set_cell([a * v * n for v, n in zip(cell, size)], scale_atoms=True)

    if vacuum is not None:
        slab.center(vacuum=vacuum, axis=2)
    
    slab.adsorbate_info['cell'] = surface_cell
    slab.adsorbate_info['sites'] = sites
    
    return slab


def fcc211(symbol, size, a=None, vacuum=None, orthogonal=True):
    """FCC(211) surface.

    Does not currently support special adsorption sites.

    Currently only implemented for *orthogonal=True* with size specified
    as (i, j, k), where i, j, and k are number of atoms in each direction.
    i must be divisible by 3 to accomodate the step width.
    """
    if not orthogonal:
        raise NotImplementedError('Only implemented for orthogonal '
                                  'unit cells.')
    if size[0] % 3 != 0:
        raise NotImplementedError('First dimension of size must be '
                                  'divisible by 3.')
    atoms = FaceCenteredCubic(symbol,
                              directions=[[1, -1, -1],
                                          [0, 2, -2],
                                          [2, 1, 1]],
                              miller=(None, None, (2, 1, 1)),
                              latticeconstant=a,
                              size=(1, 1, 1),
                              pbc=True)
    z = (size[2] + 1) // 2
    atoms = atoms.repeat((size[0] // 3, size[1], z))
    if size[2] % 2:  # Odd: remove bottom layer and shrink cell.
        remove_list = [atom.index for atom in atoms
                       if atom.z < atoms[1].z]
        del atoms[remove_list]
        dz = atoms[0].z
        atoms.translate((0., 0., -dz))
        atoms.cell[2][2] -= dz
    atoms.center(vacuum=vacuum, axis=2)
    # Renumber systematically from top down.
    orders = [(atom.index, round(atom.x, 3), round(atom.y, 3),
               -round(atom.z, 3), atom.index) for atom in atoms]
    orders.sort(key=itemgetter(3, 1, 2))
    newatoms = atoms.copy()
    for index, order in enumerate(orders):
        newatoms[index].position = atoms[order[0]].position.copy()
    return newatoms


def fcc111_root(symbol, root, size, a=None, vacuum=0.0,
                orthogonal=False, search_zone=(20, 20)):
    """FCC(111) surface maniupulated to repeat with *root*
    number of atoms in each x/y plane.

    The cell is generated as a rotated 60-120-60-120
    cell in the x/y plane.  The rotation and size is chosen
    to allow for the length of the cell vectors to be equal
    to the root of *root* with a lattice constant of 2**0.5.

    *root* should be given as a positive whole number."""
    atomic_number = atomic_numbers[symbol]
    if orthogonal:
        raise NotImplementedError('Only implemented for orthogonal '
                                  'unit cells.')
    if a is None:
        if reference_states[atomic_number]['symmetry'] == 'fcc':
            a = reference_states[atomic_number]['a']
        else:
            raise ValueError("Can't guess lattice constant for %s-%s!"
                             % ('fcc', symbol))

    searchx, searchy = search_zone

    c = 0.5
    s = (3**0.5) / 2.

    def rhomb(x, y):
        return float(x + (c * y)), float(s * y)

    desired = root**0.5
    location = None

    locations = []
    for iy in range(searchy):
        for ix in range(searchx):
            x, y = rhomb(ix, iy)
            locations.append([x, y])
    distances = np.sqrt((np.array(locations)**2).sum(1))
    for index, dist in enumerate(distances):
        if abs(dist - desired) <= 1e-13:
            location = locations[index]

    if location is None:
        raise ValueError(
            "Can't find a root cell in the searched zone of size \
             (%d, %d). A larger zone may be needed for large root values" %
            search_zone)

    angle = acos(location[0] / desired)

    cutting_board = fcc111(symbol, (searchx, searchy, size[2]), a=2**0.5)
    cutting_board.translate((-searchx, 0, 0))
    cutting_board += fcc111(symbol, (searchx, searchy, size[2]), a=2**0.5)

    cutting_board.rotate('z', -angle)

    cutting_board.set_cell(((desired, 0, 0),
                            (desired / 2, desired * ((3**0.5) / 2), 0),
                            (0, 0, 2**0.5 * size[2])), scale_atoms=False)

    cell = cutting_board.get_cell()

    remove = []
    for index, position in enumerate(cutting_board.positions):
        if not (0 < position[0] < (cell[0][0] + cell[1][0])):
            remove.append(index)
        if not (0 < position[1] < cell[1][1]):
            remove.append(index)
    del cutting_board[remove]

    def remove_doubles():
        scaled = cutting_board.get_scaled_positions()

        remove = []
        for index, position in enumerate(scaled):
            for inner_index in range(index):
                inner_position = scaled[inner_index]
                if hypot(position[0] - inner_position[0],
                         position[1] - inner_position[1]) < 1e-10:
                    remove.append(index)

        del cutting_board[remove]

        cutting_board.set_scaled_positions(
            cutting_board.get_scaled_positions())

    remove_doubles()
    cutting_board.translate((0.1, 0.1, 0))
    remove_doubles()

    cell *= a / (2 ** (0.5))
    cutting_board.set_cell(cell, scale_atoms=True)
    cutting_board.adsorbate_info = {}

    cutting_board *= (size[0], size[1], 1)
    cutting_board.center(axis=2, vacuum=vacuum)

    return cutting_board
    
    
def mx2(formula='MoS2', kind='2H', a=3.18, thickness=3.19,
        size=(1, 1, 1), vacuum=7.5):
    """Create three-layer 2D materials with hexagonal structure.
    
    For metal dichalcogenites, ect.
    
    The kind argument accepts '2H', which gives a mirror plane symmetry
    and '1T', which gives an inversion symmetry."""
    
    if kind == '2H':
        basis = [(0, 0, 0),
                 (2 / 3, 1 / 3, 0.5 * thickness),
                 (2 / 3, 1 / 3, -0.5 * thickness)]
    elif kind == '1T':
        basis = [(0, 0, 0),
                 (2 / 3, 1 / 3, 0.5 * thickness),
                 (1 / 3, 2 / 3, -0.5 * thickness)]
    else:
        raise ValueError('Structure not recognized')
    
    cell = [[a, 0, 0], [-a / 2, a * 3**0.5 / 2, 0], [0, 0, 1]]
    
    atoms = Atoms(formula, cell=cell, scaled_positions=basis, pbc=(1, 1, 0))
    atoms = atoms.repeat(size)
    atoms.center(vacuum=vacuum, axis=2)
    
    return atoms
