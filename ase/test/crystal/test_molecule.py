from ase.optimize import BFGS
from ase.atoms import Atoms
from ase.calculators.crystal import CRYSTAL

geom = Atoms('OHH',
             positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0)])

geom.set_calculator(CRYSTAL(label='water',
                            guess=True,
                            basis='sto-3g',
                            xc='PBE',
                            otherkeys=['scfdir', 'anderson',
                                       ['maxcycles', '500'],
                                       ['toldee', '6'],
                                       ['tolinteg', '7 7 7 7 14'],
                                       ['fmixing', '90']]))

opt = BFGS(geom)
opt.run(fmax=0.05)

final_energy = geom.get_potential_energy()
assert abs(final_energy + 2047.34531091) < 1.0
