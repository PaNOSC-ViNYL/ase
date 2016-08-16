from ase.collections import g2
from ase.calculators.octopus import Octopus
from ase.optimize import QuasiNewton


# Ethanol molecule with somewhat randomized initial positions:
system = g2['CH3CH2OH']
system.rattle(stdev=0.1, seed=42)
system.center(vacuum=3.0)

calc = Octopus(label='ethanol',
               Spacing=0.25,
               BoxShape='parallelepiped')
system.set_calculator(calc)

opt = QuasiNewton(system, logfile='opt.log', trajectory='opt.traj')
opt.run(fmax=0.05)
