import sys

from ase import *
from ase import Atoms
from ase.lattice.surface import *
from math import *
from ase.optimize.bfgs_dev2 import BFGS_dev
from ase.optimize.bfgs import BFGS

opt = sys.argv[1]
maxstep = float(sys.argv[2])
alpha = float(sys.argv[3])
if len(sys.argv) > 4:
    memory = int(sys.argv[4])
else:
    memory = 0
print '%-15s %5s' % ('Optimizer', '#1')
#for opt in ['LBFGS', 'BFGS', 'BFGS_dev']:
if 1:
    if opt == 'BFGS':
        QuasiNewton = BFGS
    #elif opt == 'oldBFGS':
    #    QuasiNewton = BFGS
    #elif opt == 'newLBFGS':
    #    QuasiNewton = newLBFGS
    elif opt == 'LBFGS':
        QuasiNewton = LBFGS
    #elif opt == 'BFGS_dev':
    #    QuasiNewton = BFGS_dev

    maxstep = None
    if maxstep is not None:
        name = opt + '-%.2f' % maxstep
    else:
        name = opt
 
    #Set up the adsorbate
    zpos = cos(134.3/2.0*pi/180.0)*1.197
    xpos = sin(134.3/2.0*pi/180.0)*1.19
    no2 =Atoms('CO', positions=[(-xpos+1.2,0,-zpos), (-xpos+1.2,-1.1,-zpos)])
    
    # Surface slab
    slab =fcc111('Au', size=(2, 2, 4),vacuum=2*5, orthogonal = True )
    slab.center()
    add_adsorbate(slab,no2,1.5,'bridge')
    slab.set_pbc((True,True,False))
    
    calc = EMT()
    slab.set_calculator(calc)
    #constraints
    constraint = FixAtoms(mask=[(a.tag == 4) or (a.tag == 3) or (a.tag==2) for a in slab])
    slab.set_constraint(constraint)
    dyn = QuasiNewton(slab, maxstep = 0.04, logfile = name + '-CO_Au111.log', trajectory = name + '-CO_Au111.traj', alpha = alpha)
    dyn.run(fmax=0.05, steps=200)
    print '%-15s %5i' % (opt, dyn.get_number_of_steps())
