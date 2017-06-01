import numpy as np
from ase.calculators.gulp import GULP
from ase.optimize import BFGS
from ase.build import bulk
from ase.build import molecule

#atoms = bulk('Si').repeat((4, 4, 4))
atoms = molecule('H2O')
atoms1 = atoms.copy()
atoms1.calc = GULP(library='reaxff.lib')
opt1 = BFGS(atoms1,trajectory='bfgs.traj')
opt1.run(fmax=0.005)

atoms2 = atoms.copy()
opt2 = GULP.get_optimizer(atoms2,library='reaxff.lib')
opt2.run()

print(np.abs(opt1.atoms.positions - opt2.atoms.positions))
assert np.abs(opt1.atoms.positions - opt2.atoms.positions).max() < 1e-5
