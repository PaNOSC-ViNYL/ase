from ase import Atoms
from ase.calculators.crystal import CRYSTAL
import os

a0 = 5.43
bulk = Atoms('Si2', [(0, 0, 0),
                     (0.25, 0.25, 0.25)],
             pbc=True)
b = a0 / 2
bulk.set_cell([(0, b, b),
               (b, 0, b),
               (b, b, 0)], scale_atoms=True)

bulk.set_calculator(CRYSTAL(label='Si2',
                            guess=True,
                            basis='sto-3g',
                            xc='PBE',
                            kpts=(2, 2, 2),
                            otherkeys=['SCFDIR', 'ANDERSON',
                                      ['MAXCYCLES', '500'],
                                      ['TOLDEE', '5'],
                                      ['TOLINTEG', '7 7 7 7 14'],
                                      ['FMIXING', '50']],
                            ))

final_energy = bulk.get_potential_energy()
assert abs(final_energy + 15564.787949) < 1.0

files = ['SCFOUT.LOG', 'INPUT', 'optc1',
         'FORCES.DAT', 'dffit3.dat', 'OUTPUT']

for file in files:
    try:
        os.remove(file)
    except OSError:
        pass

dir = os.getcwd()
fort = os.listdir(dir)
for file in fort:
    if file.startswith("fort."):
        os.remove(os.path.join(dir, file))
