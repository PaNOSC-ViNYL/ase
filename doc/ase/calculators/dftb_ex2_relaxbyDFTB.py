from ase.calculators.dftb import Dftb
from ase.io import write, read

from ase.build import molecule
system = molecule('H2O')
calc = Dftb(label='h2o', atoms=system,
            run_manyDftb_steps=True,
            Driver_='ConjugateGradient',
            Driver_MaxForceComponent='1E-4',
            Driver_MaxSteps=1000,
            Hamiltonian_MaxAngularMomentum_='',
            Hamiltonian_MaxAngularMomentum_O='"p"',
            Hamiltonian_MaxAngularMomentum_H='"s"')
system.set_calculator(calc)
calc.calculate(system)
final = read('geo_end.gen')
write('test.final.xyz', final)
