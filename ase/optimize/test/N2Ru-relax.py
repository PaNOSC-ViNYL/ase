import sys

from ase import Atoms, EMT
from ase.constraints import FixAtoms
from ase.optimize.bfgs import BFGS
#from ase.optimize.bfgs import oldBFGS
from ase.optimize.lbfgs import LBFGS
#from ase.optimize.newlbfgs import newLBFGS
    
opt = sys.argv[1]
maxstep = float(sys.argv[2])
alpha = float(sys.argv[3])
if len(sys.argv) > 4:
    memory = int(sys.argv[4])
else:
    memory = 0
print '%-15s %5s %5s' % ('Optimizer', '#1', '#2')
if 1:
    if opt == 'BFGS':
        QuasiNewton = BFGS
    elif opt == 'LBFGS':
        QuasiNewton = LBFGS

    maxstep = None
    if maxstep is not None:
        name = opt + '-%.2f' % maxstep
    else:
        name = opt

    a = 2.70
    c = 1.59 * a
    h = 1.85
    d = 1.10

    slab = Atoms('2Cu', [(0., 0., 0.), (1/3., 1/3., -0.5*c)], 
                 tags=(0, 1),
                 pbc=(1, 1, 0))
    slab.set_cell([(a, 0, 0),
                   (a / 2, 3**0.5 * a / 2, 0),
                   (0, 0, 1)])
    slab = slab.repeat((4, 4, 1))
    slab.set_calculator(EMT())
    mask = [a.tag == 1 for a in slab]
    slab.set_constraint(FixAtoms(mask=mask))
    dyn0 = QuasiNewton(slab, logfile = name + '-Ru_slab.log', trajectory = name + '-Ru_slab.traj', maxstep = maxstep, alpha = alpha)
    dyn0.run(fmax=0.05, steps=200)

    e_slab = slab.get_potential_energy()
    x = slab.positions[0, 2] / (c / 2) * 100

    molecule = Atoms('2N', positions=[(0., 0., h),
                                      (0., 0., h + d)])
    molecule.set_calculator(EMT())
    e_N2 = molecule.get_potential_energy()
    slab.extend(molecule)

    dyn1 = QuasiNewton(slab, logfile = name + '-N2_on_Ru.log', trajectory = name + '-N2_on_Ru.traj', maxstep = maxstep, alpha = alpha)
    dyn1.run(fmax=0.05, steps=200)
    print '%-15s %5i %5i' % (opt, dyn0.get_number_of_steps(), dyn1.get_number_of_steps())

#print 'Adsorption energy:', e_slab + e_N2 - slab.get_potential_energy()
