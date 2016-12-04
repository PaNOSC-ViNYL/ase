import numpy as np
from ase.calculators.siesta import Siesta
from ase.calculators.siesta.parameters import Specie
from ase.optimize import QuasiNewton
from ase import Atoms

atoms = Atoms('CH4', np.array([
    [0.000000, 0.000000, 0.000000],
    [0.682793, 0.682793, 0.682793],
    [-0.682793, -0.682793, 0.682790],
    [-0.682793, 0.682793, -0.682793],
    [0.682793, -0.682793, -0.682793]]))

siesta = Siesta(
    species=[
        Specie(symbol='H', excess_charge=0.1)])

atoms.set_calculator(siesta)
dyn = QuasiNewton(atoms, trajectory='h.traj')
dyn.run(fmax=0.02)
