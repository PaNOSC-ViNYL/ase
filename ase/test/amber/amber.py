"""Test that amber calculator works.

This is conditional on the existence of the $AMBERHOME/bin/sander
executable.
"""
import subprocess

from ase import Atoms
from ase.calculators.amber import Amber
from ase.test import require


require('amber')

with open('2h2o.pdb', 'w') as outfile:
    outfile.write("""\
ATOM      1  O   HOH     1       0.000   0.000   0.000  1.00  0.00           O
ATOM      3  H1  HOH     1      -0.241  -0.971   0.000  1.00  0.00           H
ATOM      6  H2  HOH     1       1.050   0.100   0.100  1.00  0.00           H
TER
ATOM      2  O   HOH     2       2.200   0.000   0.000  1.00  0.00           O
ATOM      4  H2  HOH     2       2.459   0.966   0.000  1.00  0.00           H
ATOM      5  H1  HOH     2       2.459  -0.483   0.837  1.00  0.00           H
END
""")

with open('mm.in', 'w') as outfile:
    outfile.write("""\
zero step md to get energy and force
&cntrl
imin=0, nstlim=0,  ntx=1 !0 step md
cut=100, ntb=0,          !non-periodic
ntpr=1,ntwf=1,ntwe=1,ntwx=1 ! (output frequencies)
&end
END
""")

with open('tleap.in', 'w') as outfile:
    outfile.write("""\
source leaprc.protein.ff14SB
source leaprc.gaff
source leaprc.water.tip3p
mol = loadpdb 2h2o.pdb
saveamberparm mol 2h2o.top h2o.inpcrd
quit
""")

subprocess.call('tleap -f tleap.in'.split())

atoms = Atoms('OH2OH2',
              [[-0.956, -0.121, 0],
               [-1.308, 0.770, 0],
               [0.000, 0.000, 0],
               [3.903, 0.000, 0],
               [4.215, -0.497, -0.759],
               [4.215, -0.497, 0.759]])

calc = Amber(amber_exe='sander -O ',
             infile='mm.in',
             outfile='mm.out',
             topologyfile='2h2o.top',
             incoordfile='mm.crd')
calc.write_coordinates(atoms, 'mm.crd')
atoms.set_calculator(calc)

e = atoms.get_potential_energy()
assert abs(e + 0.046799672) < 5e-3
