"""Test that amber calculator works."""
from ase.calculators.amber import Amber
import ase.io as io

atoms = io.read('2h2o.pdb')
calc = Amber(amber_exe='sander -O ',
             infile='mm.in',
             outfile='mm.out',
             topologyfile='2h2o.top',
             incoordfile='mm.crd')
calc.write_coordinates(atoms, 'mm.crd')
atoms.set_calculator(calc)
e = atoms.get_potential_energy()
assert abs(e + 0.046799672) < 5e-3
