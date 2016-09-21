"""
Run Amber test that amber calculator works. This
is conditional on the existence of the $AMBERHOME/bin/sander
executable.
"""

from ase.test.amber import installed
assert installed()


from ase.calculators.amber import Amber
from ase.optimize import BFGS
import ase.io as io

atoms = io.read('2h2o.pdb')
calc = Amber(amber_exe='sander -O ',
             infile = 'mm.in', outfile = 'mm.out',
             topologyfile = '2h2o.top', incoordfile='mm.crd')
calc.write_coordinates(atoms, 'mm.crd')
atoms.set_calculator(calc)
e = atoms.get_potential_energy()
assert abs(e + 0.046799672) < 5e-3
