import numpy as np
from ase.calculators.gulp import GULP
from ase.optimize import BFGS
from ase.build import bulk

atoms = bulk('Si').repeat((4, 4, 4))

atoms1 = atoms.copy()
atoms1.calc = GULP()
opt1 = BFGS(atoms1)
opt1.run(fmax=0.05)

atoms2 = atoms.copy()
opt2 = GULP.get_optimizer(atoms2)
opt2.run(fmax=0.05)

assert np.abs(opt1.positions - opt2.positions).max() < 1e-5
