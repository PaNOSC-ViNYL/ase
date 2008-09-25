from math import sqrt

import numpy as npy


class NeighborList:
    """Neighbor list object.

    :Parameters:
      cutoffs: list of float
          List of cutoff radii - one for each atom.
      skin: float
          If no atom has moved more than the skin-distance since the
          last call to the ``update()`` method, then the neighbor list
          can be reused.  This will save some expensive rebuilds of
          the list, but extra neighbors outside the cutoff will be
          returned.

    Example::

      nl = NeighborList([2.3, 1.7])
      nl.update()
      indices, offsets = nl.get_neighbors(0)
      
    """
    
    def __init__(self, cutoffs, skin=0.3, sorted=False):
        self.cutoffs = npy.asarray(cutoffs) + skin
        self.skin = skin
        self.sorted = sorted
        self.nupdates = 0

    def update(self, atoms):
        """Make sure the list is up to date."""
        if self.nupdates == 0:
            self.build(atoms)
            return True
        
        assert (self.cell == atoms.get_pbc()).all()

        if ((atoms.positions - self.positions)**2).sum(1).max() > self.skin**2:
            self.build(atoms)
            return True

        return False
    
    def build(self, atoms):
        """Build the list."""
        self.positions = atoms.get_positions()
        pbc = atoms.get_pbc()
        self.cell = atoms.get_cell()
        rcmax = self.cutoffs.max()
        
        icell = npy.linalg.inv(self.cell)
        scaled = npy.dot(self.positions, icell)
        scaled0 = scaled.copy()

        N = []
        for i in range(3):
            if pbc[i]:
                scaled0[:, i] %= 1.0
                v = icell[:, i]
                h = 1 / sqrt(npy.dot(v, v))
                n =  int(2 * rcmax / h) + 1
            else:
                n = 0
            N.append(n)
            
        offsets = npy.empty((len(atoms), 3), int)
        (scaled0 - scaled).round(out=offsets)
        positions0 = npy.dot(scaled0, self.cell)
        natoms = len(atoms)
        indices = npy.arange(natoms)

        self.neighbors = [npy.empty(0, int) for a in range(natoms)]
        self.displacements = [npy.empty((0, 3), int) for a in range(natoms)]
        for n1 in range(0, N[0] + 1):
            for n2 in range(-N[1], N[1] + 1):
                for n3 in range(-N[2], N[2] + 1):
                    if n1 == 0 and (n2 < 0 or n2 == 0 and n3 < 0):
                        continue
                    displacement = npy.dot((n1, n2, n3), self.cell)
                    for a in range(natoms):
                        d = positions0 + displacement - positions0[a]
                        i = indices[(d**2).sum(1) <
                                    (self.cutoffs + self.cutoffs[a])**2]
                        if n1 == 0 and n2 == 0 and n3 == 0:
                            i = i[i >= a]
                        self.neighbors[a] = npy.concatenate(
                            (self.neighbors[a], i))
                        disp = npy.empty((len(i), 3), int)
                        disp[:] = (n1, n2, n3)
                        disp += offsets[i] - offsets[a]
                        self.displacements[a] = npy.concatenate(
                            (self.displacements[a], disp))

        if self.sorted:
            for a, i in enumerate(self.neighbors):
                j = i[i < a]
                offsets = self.displacements[a][i < a]
                for b, offset in zip(j, offsets):
                    self.neighbors[b] = npy.concatenate(
                        (self.neighbors[b], [a]))
                    self.displacements[b] = npy.concatenate(
                            (self.displacements[b], [-offset]))
                self.neighbors[a] = self.neighbors[a][i >= a]
                self.displacements[a] = self.displacements[a][i >= a]
                
        self.nupdates += 1

    def get_neighbors(self, a):
        """Return neighbors of atom number a.

        A list of indices and offsets to neighboring atoms is
        returned.  The positions of the neighbor atoms can be
        calculated like this::

          indices, offsets = nl.get_neighbors(42)
          for i, offset in zip(indices, offsets):
              print atoms.positions[i] + dot(offset, atoms.get_cell())

        Notice that if get_neighbors(a) gives atom b as a neighbor,
        then get_neighbors(b) will not return a as a neighbor!"""
        
        return self.neighbors[a], self.displacements[a]
