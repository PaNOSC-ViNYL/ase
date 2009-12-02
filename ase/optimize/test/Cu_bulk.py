import sys
from ase.lattice.cubic import FaceCenteredCubic
from ase import *
from ase.optimize.bfgs import BFGS
#from ase.optimize.bfgs_dev2 import BFGS_dev

opt = sys.argv[1]
maxstep = float(sys.argv[2])
alpha = float(sys.argv[3])
if len(sys.argv) > 4:
    memory = int(sys.argv[4])
else:
    memory = 0
print '%-15s %5s' % ('Optimizer', '#1')
if opt == 'BFGS':
    QuasiNewton = BFGS
elif opt == 'LBFGS':
    QuasiNewton = LBFGS

maxstep = None
if maxstep is not None:
    name = opt + '-%.2f' % maxstep
else:
    name = opt

atoms = FaceCenteredCubic(directions=[[1,-1,0], [1,1,0], [0,0,1]],
                          size=(3,3,3), symbol='Cu', pbc=(1,1,1))
atoms.rattle(stdev=0.1,seed=42)
atoms.set_calculator(EMT())

#constraint = FixAtoms(indices=[0])
#atoms.set_constraint(constraint)
dyn = QuasiNewton(atoms,maxstep=0.04,trajectory=name + '-Cu_bulk.traj', alpha = alpha)
dyn.run(fmax=0.02)
print '%-15s %5i' % (opt, dyn.get_number_of_steps())
