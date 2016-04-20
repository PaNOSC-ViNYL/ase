from ase.lattice.surface import fcc111
from ase.constraints import FixAtoms, FixBondLengths, FixInternals

slab = fcc111("Pt", (4, 4, 4))

C1 = FixAtoms([0, 2, 4])
C2 = FixBondLengths([[0, 1], [0, 2]]) 
C3 = FixInternals(bonds = [[1, [7, 8]], [1, [8, 9]]])

slab.set_constraint([C1, C2, C3])
print slab.get_constrained_indices()
