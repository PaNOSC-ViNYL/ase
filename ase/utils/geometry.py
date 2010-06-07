# Copyright (C) 2010, Jesper Friis
# (see accompanying license files for details).

"""Utility tools for convenient creation of slabs and interfaces of
different orientations."""

import numpy as np



def gcd(seq):
    """Returns greatest common divisor of integers in *seq*."""
    def _gcd(m, n):
        while n:
            m, n = n, m%n
        return m       
    return reduce(_gcd, seq)




def get_layers(atoms, miller, tolerance=0.001):
    """Returns two arrays describing which layer each atom belongs
    to and the distance between the layers and origo. 

    Parameters:

    miller: 3 integers
        The Miller indices of the planes. Actually, any direction
        in reciprocal space works, so if a and b are two float
        vectors spanning an atomic plane, you can get all layers
        parallel to this with miller=np.cross(a,b).
    tolerance: float
        The maximum distance in Angstrom along the plane normal for
        counting two atoms as belonging to the same plane.

    Returns:

    tags: array of integres
        Array of layer indices for each atom.
    levels: array of floats
        Array of distances in Angstrom from each layer to origo.

    Example:

    >>> import numpy as np
    >>> from ase.lattice.spacegroup import crystal
    >>> atoms = crystal('Al', [(0,0,0)], spacegroup=225, cellpar=4.05)
    >>> np.round(atoms.positions, decimals=5)
    array([[ 0.   ,  0.   ,  0.   ],
           [ 0.   ,  2.025,  2.025],
           [ 2.025,  0.   ,  2.025],
           [ 2.025,  2.025,  0.   ]])
    >>> get_layers(atoms, (0,0,1))
    (array([0, 1, 1, 0]), array([ 0.   ,  2.025]))
    """
    miller = np.asarray(miller)

    metric = np.dot(atoms.cell, atoms.cell.T)
    c = np.linalg.solve(metric.T, miller.T).T
    miller_norm = np.sqrt(np.dot(c, miller))
    d = np.dot(atoms.get_scaled_positions(), miller)/miller_norm

    keys = np.argsort(d)
    ikeys = np.argsort(keys)
    mask = np.concatenate(([True], np.diff(d[keys]) > tolerance))
    tags = np.cumsum(mask)[ikeys]
    if tags.min() == 1:
        tags -= 1

    levels = d[keys][mask]
    return tags, levels




def cut(atoms, a=(1,0,0), b=(0,1,0), c=(0,0,1), origo=(0,0,0), 
        nlayers=None, extend=1.0, tolerance=0.001):
    """Cuts out a cell defined by *a*, *b*, *c* and *origo* from a
    sufficiently repeated copy of *atoms*.

    Typically, this function is used to create slabs of different
    sizes and orientations. The vectors *a*, *b* and *c* are in scaled
    coordinates and defines the returned cell and should normally be
    integer-valued in order to end up with a periodic
    structure. However, for systems with sub-translations, like fcc,
    integer multiples of 1/2 or 1/3 might also make sence for some
    directions (and will be treated correctly).

    Parameters:
    
    atoms: Atoms instance
        This should correspond to a repeatable unit cell.
    a: int | 3 floats
        The a-vector in scaled coordinates of the cell to cut out. If
        integer, the a-vector will be the scaled vector from *origo* to the
        atom with index *a*.
    b: int | 3 floats
        The b-vector in scaled coordinates of the cell to cut out. If
        integer, the b-vector will be the scaled vector from *origo* to the
        atom with index *b*.
    c: int | 3 floats
        The c-vector in scaled coordinates of the cell to cut out. If
        integer, the c-vector will be the scaled vector from *origo* to the
        atom with index *c*. Not used if *nlayers* is given.
    origo: int | 3 floats
        Position of origo of the new cell in scaled coordinates. If
        integer, the position of the atom with index *origo* is used.
    nlayers: int
        If *nlayers* is not *None*, the returned cell will have
        *nlayers* atomic layers in the c-direction. The direction of
        the c-vector will be along cross(a, b) converted to real
        space, i.e. normal to the plane spanned by a and b in
        orthorombic systems.
    extend: 1 or 3 floats
        The *extend* argument scales the effective cell in which atoms
        will be included. It must either be three floats or a single
        float scaling all 3 directions.  By setting to a value just
        above one, e.g. 1.05, it is possible to all the corner and
        edge atoms in the returned cell.  This will of cause make the
        returned cell non-repeatable, but is very usefull for
        visualisation.
    tolerance: float
        Determines what is defined as a plane.  All atoms within
        *tolerance* Angstroms from a given plane will be considered to
        belong to that plane.

    Example:

    >>> import ase
    >>> from ase.lattice.spacegroup import crystal
    >>>
    # Create an aluminium (111) slab with three layers
    #
    # First an unit cell of Al
    >>> a = 4.05
    >>> aluminium = crystal('Al', [(0,0,0)], spacegroup=225,
    ...                     cellpar=[a, a, a, 90, 90, 90])
    >>>
    # Then cut out the slab
    >>> al111 = cut(aluminium, (1,-1,0), (0,1,-1), nlayers=3)
    >>>
    # Visualisation of the skutterudite unit cell 
    #
    # Again, create a skutterudite unit cell
    >>> a = 9.04
    >>> skutterudite = crystal(
    ...     ('Co', 'Sb'), 
    ...     basis=[(0.25,0.25,0.25), (0.0, 0.335, 0.158)], 
    ...     spacegroup=204, 
    ...     cellpar=[a, a, a, 90, 90, 90])
    >>>
    # Then use *origo* to put 'Co' at the corners and *extend* to
    # include all corner and edge atoms.
    >>> s = cut(skutterudite, origo=(0.25, 0.25, 0.25), extend=1.01)
    >>> ase.view(s)  # doctest: +SKIP
    """
    atoms = atoms.copy()
    cell = atoms.cell

    if isinstance(origo, int):
        origo = atoms.get_scaled_positions()[origo]
    scaled = (atoms.get_scaled_positions() - origo)%1.0
    scaled %= 1.0 # needed to ensure that all numbers are *less* than one
    atoms.set_scaled_positions(scaled)

    if isinstance(a, int):
        a = scaled[a] - origo
    if isinstance(b, int):
        b = scaled[b] - origo
    if isinstance(c, int):
        c = scaled[c] - origo

    a = np.array(a, dtype=float)
    b = np.array(b, dtype=float)
    origo = np.array(origo, dtype=float)

    if nlayers:
        miller = np.cross(a, b) # surface normal
        # The factor 36 = 2*2*3*3 is because the elements of a and b
        # might be multiples of 1/2 or 1/3 because of lattice
        # subtranslations
        if np.all(36*miller - np.rint(36*miller)) < 1e-5:
            miller = np.rint(36*miller)
            miller /= gcd(miller)
        tags, layers = get_layers(atoms, miller, tolerance)
        while tags.max() < nlayers:
            atoms = atoms.repeat(2)
            tags, layers = get_layers(atoms, miller, tolerance)
        # Convert surface normal in reciprocal space to direction in
        # real space
        metric = np.dot(cell, cell.T)
        c = np.linalg.solve(metric.T, miller.T).T
        c *= layers[nlayers]/np.sqrt(np.dot(c, miller))

    newcell = np.dot(cell, np.array([a, b, c]))
    
    # Create a new atoms object, repeated and translated such that
    # it completely covers the new cell
    scorners_newcell = np.array([[0., 0., 0.], [0., 0., 1.], 
                                 [0., 1., 0.], [0., 1., 1.], 
                                 [1., 0., 0.], [1., 0., 1.], 
                                 [1., 1., 0.], [1., 1., 1.]])
    corners = np.dot(scorners_newcell, newcell*extend)
    scorners = np.linalg.solve(cell.T, corners.T).T
    rep = np.ceil(scorners.ptp(axis=0)).astype('int') + 1
    trans = np.dot(np.floor(scorners.min(axis=0)), cell)
    atoms = atoms.repeat(rep)
    atoms.translate(trans)
    atoms.set_cell(newcell)

    # Mask out atoms outside new cell
    stol = tolerance  # scaled tolerance, XXX
    maskcell = atoms.cell*extend
    sp = np.linalg.solve(maskcell.T, (atoms.positions).T).T
    mask = np.all(np.logical_and(-stol <= sp, sp < 1-stol), axis=1)
    atoms = atoms[mask]
    
    return atoms




def stack(atoms1, atoms2, axis=2, cell=None, fix=0.5,  
          maxstrain=0.5, distance=None):
    """Return a new Atoms instance with *atoms2* added to atoms1
    along the given axis. Periodicity in all directions is
    ensured.

    The size of the final cell is determined by *cell*, except
    that the length alongh *axis* will be the sum of
    *atoms1.cell[axis]* and *atoms2.cell[axis]*. If *cell* is None,
    it will be interpolated between *atoms1* and *atoms2*, where
    *fix* determines their relative weight. Hence, if *fix* equals
    zero, the final cell will be determined purely from *atoms1* and
    if *fix* equals one, it will be determined purely from
    *atoms2*.

    An ValueError exception will be raised if the far corner of
    the unit cell of either *atoms1* or *atoms2* is displaced more
    than *maxstrain*. Setting *maxstrain* to None, disable this
    check.

    If *distance* is provided, the atomic positions in *atoms1* and
    *atoms2* as well as the cell lengths along *axis* will be
    adjusted such that the distance between the distance between
    the closest atoms in *atoms1* and *atoms2* will equal *distance*.
    This option uses scipy.optimize.fmin() and hence require scipy
    to be installed.

    Example:

    >>> import ase
    >>> from ase.lattice.spacegroup import crystal
    >>>
    # Create an Ag(110)-Si(110) interface with three atomic layers
    # on each side. 
    >>> a_ag = 4.09
    >>> ag = crystal(['Ag'], basis=[(0,0,0)], spacegroup=225, 
    ...              cellpar=[a_ag, a_ag, a_ag, 90., 90., 90.])
    >>> ag110 = cut(ag, (0, 0, 3), (-1.5, 1.5, 0), nlayers=3)
    >>>
    >>> a_si = 5.43
    >>> si = crystal(['Si'], basis=[(0,0,0)], spacegroup=227, 
    ...              cellpar=[a_si, a_si, a_si, 90., 90., 90.])
    >>> si110 = cut(si, (0, 0, 2), (-1, 1, 0), nlayers=3)
    >>>
    >>> interface = stack(ag110, si110, maxstrain=1)
    >>> ase.view(interface)  # doctest: +SKIP
    >>>
    # Once more, this time adjusted such that the distance between
    # the closest Ag and Si atoms will be 2.3 Angstrom.
    >>> interface2 = stack(ag110, si110, 
    ...                    maxstrain=1, distance=2.3)   # doctest:+ELLIPSIS
    Optimization terminated successfully.
        ...
    >>> ase.view(interface2)  # doctest: +SKIP
    """
    atoms1 = atoms1.copy()
    atoms2 = atoms2.copy()

    c1 = np.linalg.norm(atoms1.cell[axis])
    c2 = np.linalg.norm(atoms2.cell[axis])
    if cell is None:
        cell1 = atoms1.cell.copy()
        cell2 = atoms2.cell.copy()
        cell1[axis] /= c1
        cell2[axis] /= c2
        cell = cell1 + fix*(cell2 - cell1)
    cell[axis] /= np.linalg.norm(cell[axis])
    cell1 = cell.copy()
    cell2 = cell.copy()
    cell1[axis] *= c1
    cell2[axis] *= c2

    if (maxstrain and
        (((cell1 - atoms1.cell).sum(axis=0)**2).sum() > maxstrain**2 or
         ((cell2 - atoms2.cell).sum(axis=0)**2).sum() > maxstrain**2)):
        raise ValueError('Incompatible cells.')

    sp1 = np.linalg.solve(atoms1.cell.T, atoms1.positions.T).T
    sp2 = np.linalg.solve(atoms2.cell.T, atoms2.positions.T).T
    atoms1.set_cell(cell1)
    atoms2.set_cell(cell2)
    atoms1.set_scaled_positions(sp1)
    atoms2.set_scaled_positions(sp2)

    if distance is not None:
        from scipy.optimize import fmin
        def mindist(pos1, pos2):
            n1 = len(pos1)
            n2 = len(pos2)
            idx1 = np.arange(n1).repeat(n2)
            idx2 = np.tile(np.arange(n2), n1)
            return np.sqrt(((pos1[idx1] - pos2[idx2])**2).sum(axis=1).min())
        def func(x):
            t1, t2, h1, h2 = x[0:3], x[3:6], x[6], x[7]
            pos1 = atoms1.positions + t1
            pos2 = atoms2.positions + t2
            d1 = mindist(pos1, pos2 + (h1 + 1.0)*atoms1.cell[axis])
            d2 = mindist(pos2, pos1 + (h2 + 1.0)*atoms2.cell[axis])
            return (d1 - distance)**2 + (d2 - distance)**2
        atoms1.center()
        atoms2.center()
        x0 = np.zeros((8,))
        x = fmin(func, x0)
        t1, t2, h1, h2 = x[0:3], x[3:6], x[6], x[7]
        atoms1.translate(t1)
        atoms2.translate(t2)
        atoms1.cell[axis] *= 1.0 + h1
        atoms2.cell[axis] *= 1.0 + h2

    atoms2.translate(atoms1.cell[axis])
    atoms1.cell[axis] += atoms2.cell[axis]
    atoms1.extend(atoms2)
    return atoms1


#-----------------------------------------------------------------
# Self test
if __name__ == '__main__':
    import doctest
    print 'doctest: ', doctest.testmod()
