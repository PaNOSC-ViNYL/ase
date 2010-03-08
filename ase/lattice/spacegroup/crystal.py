# Copyright (C) 2010, Jesper Friis
# (see accompanying license files for details).

"""
A module for ASE for simple creation of crystalline structures from
knowledge of the space group.

"""

import numpy as np

import ase
from ase.atoms import string2symbols
from spacegroup import Spacegroup
from cell import cellpar_to_cell

__all__ = ['crystal']



def crystal(symbols, basis, spacegroup=1, setting=1, 
            cell=None, cellpar=(1,1,1,90,90,90), 
            ab_normal=(0,0,1), a_direction=None, size=(1,1,1),
            ondublicates='warn', symprec=0.001, 
            pbc=True, **kwargs):
    """Create an Atoms instance for a conventional unit cell of a
    space group.

    Parameters:

    symbols : string | sequence of strings
        Either a sequence of element symbols or a blank-separated string
        of element symbols. E.g. ('Na', 'Cl') and 'NaCl' are equivalent.
    basis : list of scaled coordinates
        Coordinates of the non-equivalent sites in units of the 
        lattice vectors.
    spacegroup : int | string | Spacegroup instance
        Space group given either as its number in International Tables
        or as its Hermann-Mauguin symbol.
    setting : 1 | 2
        Space group setting.
    cell : 3x3 matrix
        Unit cell vectors.
    cellpar : [a, b, c, alpha, beta, gamma]
        Cell parameters with angles in degree. Is not used when `cell`
        is given. 
    ab_normal : vector
        Is used to define the orientation of the unit cell relative
        to the Cartesian system when `cell` is not given. It is the
        normal vector of the plane spanned by a and b.
    a_direction : vector
        Defines the orientation of the unit cell a vector. a will be 
        parallel to the projection of `a_direction` onto the a-b plane.
    size : 3 positive integers
        How many times the conventional unit cell should be repeated
        in each direction.
    ondublicates : 'keep' | 'replace' | 'warn' | 'error'
        Action if `basis` contain symmetry-equivalent positions:
            'keep'    - ignore additional symmetry-equivalent positions
            'replace' - reolace
            'warn'    - like 'keep', but issue an UserWarning
            'error'   - raises a SpacegroupValueError
    symprec : float
        Minimum "distance" betweed two sites in scaled coordinates
        before they are counted as the same site.
    pbc : one or three bools
        Periodic boundary conditions flags.  Examples: True,
        False, 0, 1, (1, 1, 0), (True, False, False).  Default
        is True.

    Keyword arguments:

    All additional keyword arguments are passed on to the Atoms constructor. 
    Currently, probably the most useful additional keyword arguments are
    `constraint` and `calculator`.

    Examples:

    Two diamond unit cells (space group number 227)

    >>> diamond = crystal('C', [(0,0,0)], spacegroup=227, 
    ...     cellpar=[3.57, 3.57, 3.57, 90, 90, 90], size=(2,1,1))
    >>> ase.view(diamond)  # doctest: +SKIP

    A CoSb3 skutterudite unit cell containing 32 atoms

    >>> skutterudite = crystal(('Co', 'Sb'), 
    ...     basis=[(0.25,0.25,0.25), (0.0, 0.335, 0.158)], 
    ...     spacegroup=204, cellpar=[9.04, 9.04, 9.04, 90, 90, 90])
    >>> len(skutterudite)
    32
    """
    sg = Spacegroup(spacegroup, setting)
    sites, kinds = sg.equivalent_sites(basis, 
                                       ondublicates=ondublicates, 
                                       symprec=symprec)
    symbols = parse_symbols(symbols)
    symbols = [symbols[i] for i in kinds]
    if cell is None:
        cell = cellpar_to_cell(cellpar, ab_normal, a_direction)
    sites, symbols, cell = repeat(size, sites, symbols, cell)
    atoms = ase.Atoms(symbols, 
                      scaled_positions=sites, 
                      cell=cell,
                      pbc=pbc,
                      **kwargs)
    return atoms


    
def surface(symbols, basis, spacegroup=1, setting=1, 
            normal=(0,0,1), corner=0, layers=None, uvwcell=None, vacuum=0.0,
            cell=None, cellpar=(1,1,1,90,90,90), 
            ab_normal=(0,0,1), a_direction=None, size=(1,1,1),
            ondublicates='warn', symprec=0.001, 
            pbc=True, **kwargs):
    """Creates a slab of a given structure with a well-defined surface
    normal.  Atoms instances created with this function will set the
    tags property in the same way as ase.lattice.surface.
    
    Parameters:

    symbols : string | sequence of strings
        Either a sequence of element symbols or a blank-separated string
        of element symbols. E.g. ('Na', 'Cl') and 'NaCl' are equivalent.
    basis : list of scaled coordinates
        Coordinates of the non-equivalent sites in units of the 
        lattice vectors.
    spacegroup : int | string | Spacegroup instance
        Space group given either as its number in International Tables
        or as its Hermann-Mauguin symbol.
    setting : 1 | 2
        Space group setting.
    normal : 3 integers
        Miller indices of top and bottom surfaces. Not used if `uvwcell`
        is provided.
    corner : integer
        Index of atom kind in symbols that will define the slab corners.
    layers : int
        Specified the number of atomic layers with the given surface
        normal to include in the returned slab.
    uvwcell : 3x3 array
        Three vectors defining the slab. If not provided, a minimal cell
        with the given `normal` will be used.
    vacuum : float
        Thickness of appended vacuum layer in Angstrom
    cell : 3x3 array
        Unit cell vectors.
    cellpar : [a, b, c, alpha, beta, gamma]
        Cell parameters with angles in degree. Is not used when `cell`
        is given. 
    ab_normal : vector
        Is used to define the orientation of the unit cell relative
        to the Cartesian system when `cell` is not given. It is the
        normal vector of the plane spanned by a and b.
    a_direction : vector
        Defines the orientation of the unit cell a vector. a will be 
        parallel to the projection of `a_direction` onto the a-b plane.
    size : 3 positive integers
        How many times the conventional unit cell should be repeated
        in each direction. The repitition in the c direction might be
        overwritten by the `layers` argument.
    ondublicates : 'keep' | 'replace' | 'warn' | 'error'
        Action if `basis` contain symmetry-equivalent positions:
            'keep'    - ignore additional symmetry-equivalent positions
            'replace' - reolace
            'warn'    - like 'keep', but issue an UserWarning
            'error'   - raises a SpacegroupValueError
    symprec : float
        Minimum "distance" betweed two sites in scaled coordinates
        before they are counted as the same site.
    pbc : one or three bools
        Periodic boundary conditions flags.  Examples: True,
        False, 0, 1, (1, 1, 0), (True, False, False).  Default
        is True.

    Keyword arguments:

    All additional keyword arguments are passed on to the Atoms constructor. 
    Currently, probably the most useful additional keyword arguments are
    `constraint` and `calculator`.
    
    Examples:

    Aluminium slab with the (111) surface in the xy-plane

    >>> slab = surface('Al', basis=[(0,0,0)], spacegroup=225, 
    ...     cellpar=[4,4,4,90,90,90], normal=[1,1,1], ab_normal=[0,0,1],
    ...     size=(3,3,3))

    Rutile slab with 7 (110) layers and 5 Angstrom vacuum on each side

    >>> rutile = surface(['Ti', 'O'], basis=[(0,0,0), (0.3, 0.3, 0.0)], 
    ...     spacegroup=136, cellpar=[4,4,3,90,90,90], 
    ...     normal=[1,1,0], ab_normal=[0,0,1], size=(3,3,3), 
    ...     layers=7, vacuum=10)
    >>> ase.view(rutile)  # doctest: +SKIP

    """
    normal = np.asarray(normal, dtype='int')
    
    sg = Spacegroup(spacegroup, setting)

    # Calculate uvwcell if it is not given
    if uvwcell is None:
        sp, kinds = sg.equivalent_sites([basis[corner]], 
                                        ondublicates=ondublicates, 
                                        symprec=symprec)
        uvwcell = cell_from_surface(sp, normal, symprec=symprec)

    # Figure out how many times we need to repeat and translate cell
    # in order for it to completely cover uvwcell
    corners = np.array([[0., 0., 0.], [0., 0., 1.], 
                        [0., 1., 0.], [0., 1., 1.], 
                        [1., 0., 0.], [1., 0., 1.], 
                        [1., 1., 0.], [1., 1., 1.]])
    cnrs = np.dot(corners, uvwcell)
    rep = np.ceil(cnrs.ptp(axis=0)).astype('int') + 1
    trans = np.floor(cnrs.min(axis=0))

    # Convert translated position of all atoms to uvw system and mask
    # out the sites within the uvwcell
    sites, kinds = sg.equivalent_sites(basis,
                                       ondublicates=ondublicates, 
                                       symprec=symprec)
    symbols = parse_symbols(symbols)
    symbols = [symbols[i] for i in kinds]
    if cell is None:
        cell = cellpar_to_cell(cellpar, ab_normal, a_direction)
    fsites, fsymbols, fcell = repeat(rep, sites, symbols, cell, fixed=True)
    pos = np.dot(fsites + trans, cell)
    pos[np.abs(pos) < 1e-10] = 0.0
    newcell = np.dot(uvwcell, cell)  # uvwcell in Cartesian coords.
    uvw = np.linalg.solve(newcell.T, pos.T).T
    mask = np.all(np.logical_and(0 <= uvw, uvw < 1), axis=1)

    symbols = (np.array(fsymbols)[mask]).tolist()
    positions = pos[mask]
        
    # Create Atoms instance
    slab = ase.Atoms(symbols, 
                      positions=positions,
                      cell=newcell,
                      pbc=pbc,
                      **kwargs)

    # Repeat the slab
    slab = slab.repeat(size)

    # Set the number of layers
    if layers is not None:
        tags = layertags(slab)
        if layers > tags[-1]:
            slab = slab.repeat((1, 1, (layers - 1)/tags[-1] + 1))
            tags = layertags(slab)
        if layers < tags[-1]:
            mask = [i for i,tag in enumerate(tags) if tag <= layers]
            slab = slab[mask]

    # Set tags
    slab.set_tags(layertags(slab))

    # Add vacuum
    if vacuum is not None:
        slab.center(vacuum=vacuum, axis=2)

    #slab.positions[np.abs(slab.positions) < 1e-7] = 0.0

    return slab





#-----------------------------------------------------------------------
# Help functions
#-----------------------------------------------------------------------


def repeat(size, scaled_positions, symbols=None, cell=None, fixed=False):
    """Repeat the cell `size` number of times and return the tuple
    (scaled_positions, symbols, cell) with updated values. `size`
    should be three integers. If `fixed` is true, the cell size will
    not be scaled."""
    sites = np.asarray(scaled_positions, dtype='float')

    nrep = size[0]*size[1]*size[2]
    if nrep < 1:
        raise ValueError, 'Size must be positive integers'

    if nrep > 1:
        if symbols is not None:
            symbols *= nrep

        newsites = []
        for i in xrange(size[0]):
            for j in xrange(size[1]):
                for k in xrange(size[2]):
                    newsites.append(sites + (i, j, k))
        sites = np.vstack(newsites)
        if not fixed:
            sites /= np.array([size[0], size[1], size[2]], dtype='float')
            
        if cell is not None and not fixed:
            cell = (np.asarray(cell).T*size).T

    return sites, symbols, cell



def cell_from_surface(scaled_equivalent_sites, normal, symprec=0.001):
    """Return a new minimal unit cell [u,v,w] with a given uv-surface normal.

    scaled_equivalent_sites : N x 3 float array
        Symmetry-equivalent sites in scaled coordinates, defining
        possible corners of the returned unit cell.
    normal : 3 integers
        Miller indices of the uv-surface normal.
    symprec : float
        Kind of precision for angular discrepancy between two vectors
        before they are considered parallel or perpendicular.

    Returns:
    
    [u,v,w] : 3x3 array of floats
         New unit cell scaled coordinates.
    """
    sp = np.asarray(scaled_equivalent_sites, dtype='float')
    symprec2 = symprec**2
    rep = np.array([1, 1, 1])

    while True:
        rep += 1
        newsp, newkinds, newcell = repeat(rep, sp)
        icenter = ((newsp - 0.5)**2).sum(axis=1).argmin()
        center = newsp[icenter]
        sp_other = np.vstack((newsp[:icenter], newsp[icenter+1:]))

        # get positions on and outside surface 
        mask_other = np.dot(sp_other - center, normal)**2 < symprec2
        sp_surf = sp_other[mask_other]
        sp_nos = sp_other[~mask_other]

        # Continue looping until the cell is big enough to contain two
        # atoms in addition to center on surface and one atom outside
        # the surface
        if len(sp_surf) < 2 or len(sp_nos) < 1:
            continue

        # Calculate new unit cell defined by vectors u, v and w, where
        # u and v are the two shortest vectors on the surface and w is the 
        # shortest vector not on surface
        iclosest1 = ((sp_surf - center)**2).sum(axis=1).argmin()
        closest1 = sp_surf[iclosest1]
        mask = np.sum(np.cross(closest1 - center, sp_surf - center)**2, 
                      axis=1) >= symprec2
        sp_surf2 = sp_surf[mask]
        if len(sp_surf2) < 1:
            continue
        iclosest2 = ((sp_surf2 - center)**2).sum(axis=1).argmin()
        closest2 = sp_surf2[iclosest2]
        iclosest3 = ((sp_nos - center)**2).sum(axis=1).argmin()
        closest3 = sp_nos[iclosest3]
        u = closest1 - center
        v = closest2 - center
        w = closest3 - center
        uvw = np.array([u, v, w])*rep[0]
        if np.linalg.det(uvw) < 0:
            #uvw = -uvw
            uvw[2] *= -1
        break
    return uvw



def unique_coords(coords, symprec=0.001):
    """Returns a sorted array containing the unique value in `coords`,
    where two values are considered equal if they max differ with
    `symprec`.

    Example
    >>> coords = unique_coords([0, 3, 1.5, 1e-4, 2e-7, 1.499999, 3.00001])
    >>> np.round(coords, decimals=3)
    array([ 0. ,  1.5,  3. ])
    """
    s = np.sort(coords)
    mask = np.diff(s) > symprec
    u = np.concatenate((s[:1], s[1:][mask]))
    #decimals = int(np.floor(np.log10(1.0/symprec)))
    #return np.round(u, decimals=decimals)
    return u


def layertags(atoms, firstlayer=1, symprec=0.001):
    """Returns an array of layer numbers. Atoms in the first layer
    will get index `firstlayer`. E.g. for setting the tags attribute
    in the same manner as done in lattice.surface, one can do

    >>> atoms = crystal(['Si', 'O'], basis=[(0,0,0), (0.3, 0.3, 0.0)],
    ...     spacegroup=136, cellpar=[4,4,3,90,90,90], size=(4,4,4))
    >>> atoms.set_tags(layertags(atoms))
    """
    sp = atoms.get_scaled_positions()
    levels = unique_coords(sp[:,2], symprec)
    tags = np.zeros((len(atoms),), dtype='int')
    for i, level in enumerate(levels):
        mask = np.abs(sp[:,2] - level) < symprec
        tags[mask] = i + 1
    return tags


def parse_symbols(symbols):
    """Return `sumbols` as a sequence of element symbols."""
    if isinstance(symbols, basestring):
        symbols = string2symbols(symbols)
    return symbols






#-----------------------------------------------------------------
# Self test
if __name__ == '__main__':
    import doctest
    print 'doctest: ', doctest.testmod()
