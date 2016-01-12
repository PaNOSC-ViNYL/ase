from ase.optimize import BFGS
from ase.io import read, write
from ase.calculators.emt import EMT
from ase.ga.relax_attaches import VariansBreak
import sys


fname = sys.argv[1]

print('Now relaxing {0}'.format(fname))
a = read(fname)

a.set_calculator(EMT())
dyn = BFGS(a, trajectory=None, logfile=None)
vb = VariansBreak(a, dyn)
dyn.attach(vb.write)
dyn.run(fmax=0.05)

a.info['key_value_pairs']['raw_score'] = -a.get_potential_energy()

write(fname[:-5] + '_done.traj', a)

print('Done relaxing {0}'.format(fname))
