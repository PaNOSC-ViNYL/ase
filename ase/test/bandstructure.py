from ase.build import bulk
from ase.calculators.test import FreeElectrons

a = bulk('Li')
a.calc = FreeElectrons(kpts={'path': 'GHPN', 'npoints': 200})
a.get_potential_energy()
bs = a.calc.band_structure()
bs.plot()
