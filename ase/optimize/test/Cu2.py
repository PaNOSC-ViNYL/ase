from ase import Atoms
from ase.constraints import FixAtoms

a = 2.70
c = 1.59 * a

slab = Atoms('2Cu', [(0., 0., 0.), (1/3., 1/3., -0.5*c)],
             tags=(0, 1),
             pbc=(1, 1, 0))
slab.set_cell([(a, 0, 0),
               (a / 2, 3**0.5 * a / 2, 0),
               (0, 0, 1)])
slab = slab.repeat((4, 4, 1))
mask = [a.tag == 1 for a in slab]
slab.set_constraint(FixAtoms(mask=mask))
