from ase.lattice.surface import fcc111, fcc111_root
from ase.lattice.root import root_surface, root_surface_analysis

# Make sample primitive cell
primitive = fcc111("Pt", (1, 1, 3))

# Check valid roots up to root search 5 (much higher than root 5)
valid = root_surface_analysis(primitive, 5)

# Make an easy sample to check code errors
atoms1 = root_surface(primitive, 7)

# Ensure the valid roots are the first 10 valid
# roots for this system
assert valid[:10] == [1.0, 3.0, 4.0, 7.0, 9.0,
                      12.0, 13.0, 16.0, 19.0, 21.0]

# Remake easy sample using surface function
atoms2 = fcc111_root("Pt", 7, (1, 1, 3))

assert len(atoms1) == len(atoms2)
assert (atoms1.positions == atoms2.positions).all()
assert (atoms1._cell == atoms2._cell).all()
