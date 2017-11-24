from math import sqrt

import numpy as np


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
