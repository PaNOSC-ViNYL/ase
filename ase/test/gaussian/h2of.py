from ase import Atoms
from ase.calculators.gaussian import Gaussian

atoms = Atoms('OH2F', positions=[(-1.853788, -0.071113, 0.000000),
                                 (-1.892204,  0.888768, 0.000000),
                                 (-0.888854, -0.232973, 0.000000),
                                 ( 1.765870,  0.148285, 0.000000)])

label = 'h2of-anion'
calc = Gaussian(charge=-1.0,
                basis='genecp',
                method='B3LYP',
                basisfile='@def2-tzvppd.gbs/N',
                label=label,
                output='wfx',
                extra='aim=charges',
                ioplist=['6/80=1', '6/35=4000000'],
                density='current',
                population='CHelpG, Hirshfeld',
                addsec=['%s.wfx' % label],
               )

atoms.set_calculator(calc)
atoms.get_potential_energy()
