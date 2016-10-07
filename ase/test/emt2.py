from ase.calculators.emt import EMT
from ase import Atoms
from ase.build import molecule
# The fixed_cutoff parameter is needed to ensure consistent
# cutoffs in all three calculations.
a1 = Atoms('Au', calculator=EMT(fixed_cutoff=True))
e1 = a1.get_potential_energy()
a2 = molecule('C6H6', calculator=EMT(fixed_cutoff=True))
e2 = a2.get_potential_energy()
a1.translate((0, 0, 50))
a3 = a1 + a2
a3.calc = EMT(fixed_cutoff=True)
e3 = a3.get_potential_energy()
print(e1, e2, e3, e3 - e1 - e2)
assert abs(e3 - e1 - e2) < 1e-13
