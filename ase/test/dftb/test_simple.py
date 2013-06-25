""" test run for DFTB+ calculator """

from ase.test import NotAvailable

from ase.calculators.dftb import Dftb

if Dftb().get_command() is None:
    raise NotAvailable('Dftb+ required')

import os

from ase import Atoms
from ase.optimize import QuasiNewton
from ase.io import write

from ase.data.molecules import molecule
test = molecule('H2O')
test.set_calculator(Dftb(label='h2o',atoms=test,
Hamiltonian_MaxAngularMomentum_ = '',
Hamiltonian_MaxAngularMomentum_O = '"p"',
Hamiltonian_MaxAngularMomentum_H = '"s"',
))

dyn = QuasiNewton(test, trajectory='test.traj')
dyn.run(fmax=0.01)
final_energy = test.get_potential_energy()
assert abs(final_energy + 111.141945) < 5e-3
files = ['band.out', 'detailed.out', 'dftb_in.hsd', 'dftb_pin.hsd', \
    'geo_end.gen', 'geo_end.xyz', 'h2o.out', 'results.tag', 'test.traj', \
    'test.traj.bak' ]
for file in files:
    try:
        os.remove(file)
    except OSError:
        pass
