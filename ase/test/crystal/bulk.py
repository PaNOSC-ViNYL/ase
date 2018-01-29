from ase import Atoms
from ase.calculators.crystal import CRYSTAL

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
                            otherkeys=['scfdir', 'anderson',
                                       ['maxcycles', '500'],
                                       ['toldee', '6'],
                                       ['tolinteg', '7 7 7 7 14'],
                                       ['fmixing', '50']]))

final_energy = bulk.get_potential_energy()
assert abs(final_energy + 15564.787949) < 1.0
