from io import StringIO

from ase.io import read
from ase.calculators.emt import EMT
from ase.optimize.berny import Berny


atoms = read(StringIO("""\
13
Formula: Ag13
Ag             0.0             0.0             0.0
Ag             0.0         -2.0345         -2.0345
Ag          2.0345             0.0         -2.0345
Ag          2.0345         -2.0345             0.0
Ag         -2.0345             0.0         -2.0345
Ag             0.0          2.0345         -2.0345
Ag         -2.0345         -2.0345             0.0
Ag          2.0345          2.0345             0.0
Ag             0.0         -2.0345          2.0345
Ag          2.0345             0.0          2.0345
Ag         -2.0345          2.0345             0.0
Ag         -2.0345             0.0          2.0345
Ag             0.0          2.0345          2.0345
"""), format='xyz')


atoms.calc = EMT()
opt = Berny(atoms, dihedral=False)
opt.run(fmax=0.05, steps=200)
