from ase.units import Ry
from ase.calculators.siesta.parameters import Specie, PAOBasisBlock
from ase.calculators.siesta.siesta import Siesta
from ase.optimize import QuasiNewton
from ase import Atoms
import numpy as np

bud = Atoms('CH4', np.array([
    [0.000000, 0.000000, 0.100000],
    [0.682793, 0.682793, 0.682793],
    [-0.682793, -0.682793, 0.68279],
    [-0.682793, 0.682793, -0.682793],
    [0.682793, -0.682793, -0.682793]]),
    cell=[10, 10, 10],
)
# Uncomment to use the last image of the relaxation trajectory.
#bud = read('bud.traj')

c_basis = """2 nodes 1.00
0 1 S 0.20 P 1 0.20 6.00
5.00
1.00
1 2 S 0.20 P 1 E 0.20 6.00
6.00 5.00
1.00 0.95"""

specie = Specie(symbol='C', basis_set=PAOBasisBlock(c_basis))
calc = Siesta(
    label='ch4',
    basis_set='SZ',
    xc='LYP',
    mesh_cutoff=300 * Ry,
    species=[specie],
    DM_Tolerance=1e-5,
    DM_MixingWeight=0.15,
    DM_NumberPulay=3,
    ElectronicTemperature=(300, 'K'),
)

bud.set_calculator(calc)
dyn = QuasiNewton(bud, trajectory='bud.traj')
dyn.run(fmax=0.02)
e = bud.get_potential_energy()
