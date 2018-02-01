from math import sqrt

import numpy as np


def mic(dr, cell, pbc=None):
    """
    Apply minimum image convention to an array of distance vectors.

    Parameters
    ----------
    dr : array_like
        Array of distance vectors.
    cell : array_like
        Simulation cell.
    pbc : array_like, optional
        Periodic boundary conditions in x-, y- and z-direction. Default is to
        assume periodic boundaries in all directions.

    Returns
    -------
    dr : array
        Array of distance vectors, wrapped according to the minimum image
        convention.
    """
    # Check where distance larger than 1/2 cell. Particles have crossed
    # periodic boundaries then and need to be unwrapped.
    rec = np.linalg.inv(cell)
    if pbc is not None:
        rec *= np.array(pbc, dtype=int).reshape(3,1)
    dri = np.round(np.dot(dr, rec))

    # Unwrap
    return dr - np.dot(dri, cell)


def neighbor_list(quantities, a, cutoff):
    """
    Compute a neighbor list for an atomic configuration.

    Parameters
    ----------
    quantities : str
        Quantities to compute by the neighbor list algorithm. Each character
        in this string defines a quantity. They are returned in a tuple of
        the same order. Possible quantities are
            'i' : first atom index
            'j' : second atom index
            'd' : absolute distance
            'D' : distance vector
            'S' : shift vector (number of cell boundaries crossed by the bond
                  between atom i and j). With the shift vector S, the
                  distances d between can be computed as:
                  D = a.positions[j]-a.positions[i]+S.dot(a.cell)
    a : ase.Atoms
        Atomic configuration.
    cutoff : float or dict
        Cutoff for neighbor search. If single float is given, a global cutoff
        is used for all elements. A dictionary specifies cutoff for element
        pairs. Specification accepts element numbers of symbols.
        Example: {(1, 6): 1.1, (1, 1): 1.0, ('C', 'C'): 1.85}

    Returns
    -------
    i, j, ... : array
        Tuple with arrays for each quantity specified above. Indices in `i`
        are returned in ascending order 0..len(a), but the order of (i,j)
        pairs is not guaranteed.

    Examples
    --------
    Examples assume Atoms object *a* and numpy imported as *np*.
    1. Coordination counting:
        i = neighbor_list('i', a, 1.85)
        coord = np.bincount(i)

    2. Coordination counting with different cutoffs for each pair of species
        i = neighbor_list('i', a,
                           {('H', 'H'): 1.1, ('C', 'H'): 1.3, ('C', 'C'): 1.85})
        coord = np.bincount(i)

    3. Pair distribution function:
        d = neighbor_list('d', a, 10.00)
        h, bin_edges = np.histogram(d, bins=100)
        pdf = h/(4*np.pi/3*(bin_edges[1:]**3 - bin_edges[:-1]**3)) * a.get_volume()/len(a)

    4. Pair potential:
        i, j, d, D = neighbor_list('ijdD', a, 5.0)
        energy = (-C/d**6).sum()
        pair_forces = (6*C/d**5  * (D/d).T).T
        forces_x = np.bincount(j, weights=pair_forces[:, 0], minlength=len(a)) - \
                   np.bincount(i, weights=pair_forces[:, 0], minlength=len(a))
        forces_y = np.bincount(j, weights=pair_forces[:, 1], minlength=len(a)) - \
                   np.bincount(i, weights=pair_forces[:, 1], minlength=len(a))
        forces_z = np.bincount(j, weights=pair_forces[:, 2], minlength=len(a)) - \
                   np.bincount(i, weights=pair_forces[:, 2], minlength=len(a))

    5. Dynamical matrix for a pair potential stored in a block sparse format:
        from scipy.sparse import bsr_matrix
        i, j, dr, abs_dr = neighbor_list('ijDd', atoms)
        energy = (dr.T / abs_dr).T
        dynmat = -(dde * (energy.reshape(-1, 3, 1) * energy.reshape(-1, 1, 3)).T).T \
                 -(de / abs_dr * (np.eye(3, dtype=energy.dtype) - \
                   (energy.reshape(-1, 3, 1) * energy.reshape(-1, 1, 3))).T).T
        dynmat_bsr = bsr_matrix((dynmat, j, first_i), shape=(3*len(a), 3*len(a)))

        Ddiag_icc = np.empty((len(a), 3, 3))
        for x in range(3):
            for y in range(3):
                Ddiag_icc[:, x, y] = -np.bincount(i, weights=dynmat[:, x, y])

        dynmat_bsr += bsr_matrix((Ddiag_icc, np.arange(len(a)), np.arange(len(a) + 1)),
                                 shape=(3 * len(a), 3 * len(a)))

    """
    # Store pbc
    pbc = a.pbc

    # Compute reciprocal lattice vectors.
    b1_c, b2_c, b3_c = np.linalg.inv(a.cell).T

    # Compute distances of cell faces.
    face_dist_c = np.array([1 / np.linalg.norm(b1_c),
                            1 / np.linalg.norm(b2_c),
                            1 / np.linalg.norm(b3_c)])

    # Compute number of bins such that a sphere of radius cutoff fit into eight
    # neighboring bins.
    nbins_c = np.maximum((face_dist_c / cutoff).astype(int), [1, 1, 1])
    nbins = np.prod(nbins_c)

    # Compute over how many cell we need to loop in the neighbor list search.
    ndx, ndy, ndz = np.ceil(cutoff * nbins_c / face_dist_c).astype(int)

    # Sort atoms into bins.
    spos_ic = a.get_scaled_positions()
    bin_index_ic = (spos_ic*nbins_c).astype(int)
    for c in range(3):
        if pbc[c]:
            bin_index_ic[:, c] %= nbins_c[c]
        else:
            bin_index_ic[:, c] = np.clip(bin_index_ic[:, c], 0, nbins_c[c]-1)

    # Convert Cartesian bin index to unique scalar bin index.
    bin_index_i = bin_index_ic[:, 0] + \
                  nbins_c[0] * (bin_index_ic[:, 1] + \
                                nbins_c[1] * bin_index_ic[:, 2])

    # Sort by bin index
    i = np.argsort(bin_index_i)
    # atom_i contains atom index in new sort order.
    atom_i = np.arange(len(a))[i]
    bin_index_i = bin_index_i[i]

    # Find max number of atoms per bin
    max_nat_per_bin = np.bincount(bin_index_i).max()

    # Sort atoms into bins: bins_ba contains for each bin (identified by its
    # scalar bin index) a list of atoms inside that bin. This list is
    # homogeneous, i.e. has the same size *max_nat_per_bin* for all bins.
    # The list is padded with -1 values.
    bins_ba = -np.ones([nbins, max_nat_per_bin], dtype=int)
    for i in range(max_nat_per_bin):
        # Create a mask array that identifies the first atom of each bin.
        m = np.append([True], bin_index_i[:-1] != bin_index_i[1:])
        # Assign all first atoms.
        bins_ba[bin_index_i[m], i] = atom_i[m]

        # Remove atoms that we just sorted into bins_ba. The next "first"
        # atom will be the second and so on.
        m = np.logical_not(m)
        atom_i = atom_i[m]
        bin_index_i = bin_index_i[m]

    # Make sure that all atoms have been sorted into bins.
    assert len(atom_i) == 0
    assert len(bin_index_i) == 0

    # Now we construct neighbor pairs by pairing up all atoms within a bin or
    # between bin and neighboring bin. ix_pn is a helper buffer that contains
    # all potential pairs of atoms between two bins, i.e. it is a list of
    # length max_nat_per_bin**2.
    ix_pn = np.indices((max_nat_per_bin, max_nat_per_bin), dtype=int)
    ix_pn = ix_pn.reshape(2, -1)

    # Initialized empty neighbor list buffers.
    i_n = []
    j_n = []
    Sx_n = []
    Sy_n = []
    Sz_n = []

    # This is the main neighbor list search. We loop over neighboring bins and
    # then construct all possible pairs of atoms between two bins, assuming
    # that each bin contains exactly max_nat_per_bin atoms. We then throw
    # out pairs involving pad atoms with atom index -1 below.
    bz_xyz, by_xyz, bx_xyz = np.meshgrid(np.arange(nbins_c[2]),
                                         np.arange(nbins_c[1]),
                                         np.arange(nbins_c[0]), indexing='ij')
    # The memory layout of bx_xyz, by_xyz, bz_xyz is such that computing the
    # respective bin index leads to a linearly increasing consecutive list.
    # The following assert statement succeeds:
    #     b_b = (bx_xyz + nbins_c[0] * (by_xyz + nbins_c[1] * bz_xyz)).ravel()
    #     assert (b_b == np.arange(np.prod(nbins_c))).all()
    for dz in range(-ndz, ndz+1):
        for dy in range(-ndy, ndy+1):
            for dx in range(-ndx, ndx+1):
                # First atom in pair.
                i_n += [bins_ba[:, ix_pn[0]]]

                # Bin index of neighboring bin and shift vector.
                sx_xyz, bx1_xyz = np.divmod(bx_xyz + dx, nbins_c[0])
                sy_xyz, by1_xyz = np.divmod(by_xyz + dy, nbins_c[1])
                sz_xyz, bz1_xyz = np.divmod(bz_xyz + dz, nbins_c[2])
                b1_b = (bx1_xyz + nbins_c[0] * 
                    (by1_xyz + nbins_c[1] * bz1_xyz)).ravel()

                # Second atom in pair.
                j_n += [bins_ba[b1_b][:, ix_pn[1]]]

                # Shift vectors.
                Sx_n += [sx_xyz.reshape(-1, 1)]
                Sy_n += [sy_xyz.reshape(-1, 1)]
                Sz_n += [sz_xyz.reshape(-1, 1)]

    # Flatten overall neighbor list.
    i_n = np.ravel(i_n)
    j_n = np.ravel(j_n)
    Sx_n = np.ravel(Sx_n)
    Sy_n = np.ravel(Sy_n)
    Sz_n = np.ravel(Sz_n)
    # Spread out shift vectors overs pairs.
    Sx_n = np.resize(Sx_n, (max_nat_per_bin**2, len(Sx_n))).T.ravel()
    Sy_n = np.resize(Sy_n, (max_nat_per_bin**2, len(Sy_n))).T.ravel()
    Sz_n = np.resize(Sz_n, (max_nat_per_bin**2, len(Sz_n))).T.ravel()
    S_n = np.transpose([Sx_n, Sy_n, Sz_n])

    #
    # Now we need to remove lots of stuff from our neighbor list...
    #

    # We have created too many pairs because we assumed each bin has exactly
    # max_nat_per_bin atoms. Remove all surperfluous pairs. Those are pairs
    # that involve an atom with index -1.
    m = np.logical_and(i_n != -1, j_n != -1)
    i_n = i_n[m]
    j_n = j_n[m]
    S_n = S_n[m]

    # Remove all self-pairs that do not cross the cell boundary.
    m = np.logical_not(np.logical_and(i_n == j_n, (S_n == 0).all(axis=1)))
    i_n = i_n[m]
    j_n = j_n[m]
    S_n = S_n[m]

    # For nonperiodic directions, remove any bonds that cross the domain
    # boundary.
    for c in range(3):
        if not pbc[c]:
            m = S_n[:, c] == 0
            i_n = i_n[m]
            j_n = j_n[m]
            S_n = S_n[m]

    # Sort neighbor list.
    i = np.argsort(i_n)
    i_n = i_n[i]
    j_n = j_n[i]
    S_n = S_n[i]

    # Compute distance vectors.
    dr_nc = a.positions[j_n] - a.positions[i_n] + S_n.dot(a.cell)
    abs_dr_n = np.sqrt(np.sum(dr_nc*dr_nc, axis=1))

    # We have still created too many pairs. Only keep those with distance
    # smaller than cutoff.
    m = abs_dr_n < cutoff
    i_n = i_n[m]
    j_n = j_n[m]
    S_n = S_n[m]
    dr_nc = dr_nc[m]
    abs_dr_n = abs_dr_n[m]

    # Assemble return tuple.
    retvals = []
    for q in quantities:
        if q == 'i':
            retvals += [i_n]
        elif q == 'j':
            retvals += [j_n]
        elif q == 'D':
            retvals += [dr_nc]
        elif q == 'd':
            retvals += [abs_dr_n]
        elif q == 'S':
            retvals += [S_n]
        else:
            raise ValueError('Unsupported quantity specified.')
    if len(retvals) == 1:
        return retvals[0]
    else:
        return tuple(retvals)


def first_neighbors(nat, i):
    """
    Compute an index array pointing to the ranges within the neighbor list that
    contain the neighbors for a certain atom.

    Parameters
    ----------
    nat : int
        Total number of atom.
    i : array_like
        First atom index 'i' of the neighbor list.

    Returns
    -------
    s : array
        Array containing pointers to the start and end location of the neighbors
        of a certain atom. Neighbors of atom k have indices from s[k] to
        s[k+1]-1.
    """
    s = -np.ones(nat+1, dtype=int)
    i = np.asarray(i)
    m = i[:-1] != i[1:]
    s[i[0]] = 0
    s[-1] = len(i)
    s[i[1:][m]] = (np.arange(len(m))+1)[m]
    m = s == -1
    while m.any():
        s[m] = s[np.arange(nat+1)[m]+1]
        m = s == -1
    return s


class NeighborList:
    """Neighbor list object.

    cutoffs: list of float
        List of cutoff radii - one for each atom. If the spheres (defined by
        their cutoff radii) of two atoms overlap, they will be counted as
        neighbors.
    skin: float
        If no atom has moved more than the skin-distance since the
        last call to the ``update()`` method, then the neighbor list
        can be reused.  This will save some expensive rebuilds of
        the list, but extra neighbors outside the cutoff will be
        returned.
    self_interaction: bool
        Should an atom return itself as a neighbor?
    bothways: bool
        Return all neighbors.  Default is to return only "half" of
        the neighbors.

    Example::

      nl = NeighborList([2.3, 1.7])
      nl.update(atoms)
      indices, offsets = nl.get_neighbors(0)
    """

    def __init__(self, cutoffs, skin=0.3, sorted=False, self_interaction=True,
                 bothways=False):
        self.nl = PrimitiveNeighborList(cutoffs, skin, sorted,
                                        self_interaction, bothways)

    def update(self, atoms):
        return self.nl.update(atoms.pbc, atoms.cell, atoms.positions)

    def get_neighbors(self, a):
        return self.nl.get_neighbors(a)

    @property
    def nupdates(self):
        return self.nl.nupdates

    @property
    def nneighbors(self):
        return self.nl.nneighbors

    @property
    def npbcneighbors(self):
        return self.nl.npbcneighbors

class PrimitiveNeighborList:
    """Neighbor list that works without Atoms objects.

    This is less fancy, but can be used to avoid conversions between
    scaled and non-scaled coordinates which may affect cell offsets
    through rounding errors.
    """
    def __init__(self, cutoffs, skin=0.3, sorted=False, self_interaction=True,
                 bothways=False, use_scaled_positions=False):
        self.cutoffs = np.asarray(cutoffs) + skin
        self.skin = skin
        self.sorted = sorted
        self.self_interaction = self_interaction
        self.bothways = bothways
        self.nupdates = 0
        self.use_scaled_positions = use_scaled_positions
        self.nneighbors = 0
        self.npbcneighbors = 0

    def update(self, pbc, cell, coordinates):
        """Make sure the list is up to date."""

        if self.nupdates == 0:
            self.build(pbc, cell, coordinates)
            return True

        if ((self.pbc != pbc).any() or
            (self.cell != cell).any() or
            ((self.coordinates - coordinates)**2).sum(1).max() >
            self.skin**2):
            self.build(pbc, cell, coordinates)
            return True

        return False

    def build(self, pbc, cell, coordinates):
        """Build the list.

        Coordinates are taken to be scaled or not according
        to self.use_scaled_positions.
        """
        self.pbc = pbc = np.array(pbc, copy=True)
        self.cell = cell = np.array(cell, copy=True)
        self.coordinates = coordinates = np.array(coordinates, copy=True)

        if len(self.cutoffs) != len(coordinates):
            raise ValueError('Wrong number of cutoff radii: {0} != {1}'
                             .format(len(self.cutoffs), len(coordinates)))

        if len(self.cutoffs) > 0:
            rcmax = self.cutoffs.max()
        else:
            rcmax = 0.0

        icell = np.linalg.pinv(cell)

        if self.use_scaled_positions:
            scaled = coordinates
            positions = np.dot(scaled, cell)
        else:
            positions = coordinates
            scaled = np.dot(positions, icell)

        scaled0 = scaled.copy()

        N = []
        for i in range(3):
            if self.pbc[i]:
                scaled0[:, i] %= 1.0
                v = icell[:, i]
                h = 1 / sqrt(np.dot(v, v))
                n = int(2 * rcmax / h) + 1
            else:
                n = 0
            N.append(n)

        offsets = (scaled0 - scaled).round().astype(int)
        positions0 = positions + np.dot(offsets, self.cell)
        natoms = len(positions)
        indices = np.arange(natoms)

        self.nneighbors = 0
        self.npbcneighbors = 0
        self.neighbors = [np.empty(0, int) for a in range(natoms)]
        self.displacements = [np.empty((0, 3), int) for a in range(natoms)]
        for n1 in range(0, N[0] + 1):
            for n2 in range(-N[1], N[1] + 1):
                for n3 in range(-N[2], N[2] + 1):
                    if n1 == 0 and (n2 < 0 or n2 == 0 and n3 < 0):
                        continue
                    displacement = np.dot((n1, n2, n3), self.cell)
                    for a in range(natoms):
                        d = positions0 + displacement - positions0[a]
                        i = indices[(d**2).sum(1) <
                                    (self.cutoffs + self.cutoffs[a])**2]
                        if n1 == 0 and n2 == 0 and n3 == 0:
                            if self.self_interaction:
                                i = i[i >= a]
                            else:
                                i = i[i > a]
                        self.nneighbors += len(i)
                        self.neighbors[a] = np.concatenate(
                            (self.neighbors[a], i))
                        disp = np.empty((len(i), 3), int)
                        disp[:] = (n1, n2, n3)
                        disp += offsets[i] - offsets[a]
                        self.npbcneighbors += disp.any(1).sum()
                        self.displacements[a] = np.concatenate(
                            (self.displacements[a], disp))

        if self.bothways:
            neighbors2 = [[] for a in range(natoms)]
            displacements2 = [[] for a in range(natoms)]
            for a in range(natoms):
                for b, disp in zip(self.neighbors[a], self.displacements[a]):
                    neighbors2[b].append(a)
                    displacements2[b].append(-disp)
            for a in range(natoms):
                nbs = np.concatenate((self.neighbors[a], neighbors2[a]))
                disp = np.array(list(self.displacements[a]) +
                                displacements2[a])
                # Force correct type and shape for case of no neighbors:
                self.neighbors[a] = nbs.astype(int)
                self.displacements[a] = disp.astype(int).reshape((-1, 3))

        if self.sorted:
            for a, i in enumerate(self.neighbors):
                mask = (i < a)
                if mask.any():
                    j = i[mask]
                    offsets = self.displacements[a][mask]
                    for b, offset in zip(j, offsets):
                        self.neighbors[b] = np.concatenate(
                            (self.neighbors[b], [a]))
                        self.displacements[b] = np.concatenate(
                            (self.displacements[b], [-offset]))
                    mask = np.logical_not(mask)
                    self.neighbors[a] = self.neighbors[a][mask]
                    self.displacements[a] = self.displacements[a][mask]

        self.nupdates += 1

    def get_neighbors(self, a):
        """Return neighbors of atom number a.

        A list of indices and offsets to neighboring atoms is
        returned.  The positions of the neighbor atoms can be
        calculated like this::

          indices, offsets = nl.get_neighbors(42)
          for i, offset in zip(indices, offsets):
              print(atoms.positions[i] + dot(offset, atoms.get_cell()))

        Notice that if get_neighbors(a) gives atom b as a neighbor,
        then get_neighbors(b) will not return a as a neighbor - unless
        bothways=True was used."""

        return self.neighbors[a], self.displacements[a]
