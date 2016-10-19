from math import sin, cos, pi
from ase import Atoms
from ase.constraints import FixAtoms
from ase.build import fcc111, add_adsorbate

zpos = cos(134.3/2.0*pi/180.0)*1.197
xpos = sin(134.3/2.0*pi/180.0)*1.19
co = Atoms('CO', positions=[(-xpos+1.2,0,-zpos), (-xpos+1.2,-1.1,-zpos)])

# Surface slab
slab =fcc111('Au', size=(2, 2, 4),vacuum=2*5, orthogonal = True )
slab.center()
add_adsorbate(slab, co,1.5,'bridge')
slab.set_pbc((True,True,False))

#constraints
constraint = FixAtoms(mask=[(a.tag == 4) or (a.tag == 3) or (a.tag==2) for a in slab])
slab.set_constraint(constraint)
del co
