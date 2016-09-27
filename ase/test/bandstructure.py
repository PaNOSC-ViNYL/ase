from ase.build import bulk
from ase.calculators.test import FreeElectrons
from ase.dft.kpoints import special_paths
from ase.dft.band_structure import BandStructure

a = bulk('Cu')
path = special_paths['fcc']
a.calc = FreeElectrons(kpts={'path': path, 'npoints': 200})
a.get_potential_energy()
bs = a.calc.band_structure()
bs.plot(emax=10, filename='bs.png', show=False)
bs.write('hmm.pckl')
bs2 = BandStructure(filename='hmm.pckl')
