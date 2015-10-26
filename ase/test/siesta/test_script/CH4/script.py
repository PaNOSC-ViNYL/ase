from ase.io import read
from ase.calculators.siesta import Siesta
from ase import Atoms
import numpy as np

bud = Atoms('CH4', np.array([
          [0.000000,  0.000000,  0.000000],
          [0.682793,  0.682793,  0.682793],
          [-0.682793, -0.682793,  0.68279],
          [-0.682793,  0.682793, -0.682793],
          [0.682793, -0.682793, -0.682793]]))

calc = Siesta(label='ch4-d')

calc.MS.set_siesta_param('SystemName', 'METHANE')
calc.MS.set_siesta_param('PAO.BasisSize', 'SZ')
calc.MS.set_siesta_param('LatticeConstant', 10, unit='Ang') 
calc.MS.set_siesta_param('AtomicCoordinatesFormat', 'Ang')
calc.MS.set_siesta_param('AtomCoorFormatOut', 'Ang')
calc.MS.set_siesta_param('AtomicCoordinatesOrigin', np.array([0.127, 0.745, -0.33]))
calc.MS.set_siesta_param('XC.Functional', 'LDA')
calc.MS.set_siesta_param('XC.Authors', 'CA')
calc.MS.set_siesta_param('MeshCutoff', 30, unit='Ry')
calc.MS.set_siesta_param('MaxSCFIterations', 50)
calc.MS.set_siesta_param('DM.MixingWeight', 0.15)
calc.MS.set_siesta_param('DM.NumberPulay', 3)
calc.MS.set_siesta_param('DM.Tolerance', 1E-4)
calc.MS.set_siesta_param('SolutionMethod', 'diagon')
calc.MS.set_siesta_param('ElectronicTemperature', 25, unit = 'meV')
calc.MS.set_siesta_param('LongOutput', '.true.')


bud.set_calculator(calc)

e = bud.get_potential_energy()
