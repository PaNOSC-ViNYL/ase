from ase.io import read
from ase.calculators.siesta import Siesta
from ase import Atoms
import numpy as np

bud = Atoms('CO2', [(0.0, 0.0, 0.0), (-1.178, 0.0, 0.0), (1.178, 0.0, 0.0)])
calc = Siesta(label='siesta')

calc.MS.set_siesta_param('MeshCutOff', 250)
calc.MS.set_siesta_param('DM.MixingWeight', 0.08)
calc.MS.set_siesta_param('DM.NumberPulay', 3)
calc.MS.set_siesta_param('DM.NumberKick', 20)
calc.MS.set_siesta_param('DM.KickMixingWeight', 0.15)
calc.MS.set_siesta_param('SolutionMethod', 'Diagon')
calc.MS.set_siesta_param('MaxSCFIterations', 500)
calc.MS.set_siesta_param("PAO.BasisSize", "DZP")
calc.MS.set_siesta_param('PAO.BasisType', 'split')
calc.MS.set_siesta_param('PAO.EnergyShift', 0.003674931, unit = 'Ry')
calc.MS.set_siesta_param('WriteCoorXmol', '.True.')
calc.MS.set_siesta_param('COOP.Write', '.True.')
calc.MS.set_siesta_param('WriteDenchar', '.True.')

calc.MS.set_siesta_param("DM.Tolerance", 0.0001)
calc.MS.set_siesta_param("DM.MixingWeight", 0.01)
calc.MS.set_siesta_param("MaxSCFIterations", 500)
calc.MS.set_siesta_param("DM.NumberPulay", 4)
calc.MS.set_siesta_param("MD.TypeOfRun", "CG")
calc.MS.set_siesta_param("MD.NumCGsteps", 0)
calc.MS.set_siesta_param("MD.MaxForceTol", 0.02, unit="eV/Ang")

bud.set_calculator(calc)

e = bud.get_potential_energy()
