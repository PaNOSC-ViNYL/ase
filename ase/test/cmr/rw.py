import numpy as np

def array_almost_equal(a1, a2, tol=np.finfo(type(1.0)).eps):
    """Replacement for old numpy.testing.utils.array_almost_equal."""
    return (np.abs(a1 - a2) < tol).all()

from ase.test import NotAvailable
# this test should be run with cmr!
try:
    import cmr
except ImportError:
    raise NotAvailable('CMR version>0.3.2 is required')

from ase.calculators.emt import EMT

from ase.io import read, write

from ase.structure import molecule

m1 = molecule('O2')
m1.center(2.0)

try:
    cmr.atoms2cmr(m1).write("O2.db")
except:
    raise NotAvailable('CMR version>0.3.2 is required')

m1.set_calculator(EMT())
e1 = m1.get_potential_energy()
f1 = m1.get_forces()

m2 = read("O2.db")

m2.set_calculator(EMT())
e2 = m2.get_potential_energy()
f2 = m1.get_forces()

# assume atoms definitions are the same if energy/forces are the same: can we do better?
assert abs(e1-e2) < 1.e-6, str(e1) + ' ' + str(e2)
assert array_almost_equal(f1, f2, tol=1.e-6)
