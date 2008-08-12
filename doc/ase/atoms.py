# creates: Au-wire.png
from ase import *
from ase.io.pov import povpng
d = 2.9
L = 10.0
wire = Atoms('Au',
             positions=[(0, L / 2, L / 2)],
             cell=(d, L, L),
             pbc=(1, 0, 0))
wire *= (6, 1, 1)
wire.positions[:, 0] -= 2 * d
wire.cell[0, 0] = d
#view(wire, block=1)
povpng('Au-wire.png', wire,
       show_unit_cell=2,
       rotation='12x,6y')
