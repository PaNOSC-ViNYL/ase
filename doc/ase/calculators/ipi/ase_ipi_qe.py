from ase.build import molecule
from ase.calculators.espresso import Espresso
from ase.units import Bohr
from ase.optimize import BFGS

atoms = molecule('H2O')
atoms.center(vacuum=3.0)

#command = 'srun pw.x < {prefix}.pwi -ipi {host}:{port} > {prefix}.pwo'

command = 'pw.x < {prefix}.pwi -ipi {host}:{port} > {prefix}.pwo'
espresso = Espresso(command=command,
                    ecutwfc=40.0,
                    pseudopotentials={'H': 'H.pbe-rrkjus.UPF',
                                      'O': 'O.pbe-rrkjus.UPF'},
                    pseudo_dir='.')

opt = BFGS(atoms, trajectory='opt.traj',
           logfile='opt.log')

with espresso.ipi(host='localhost', port=31415) as calc:
    atoms.calc = calc
    opt.run(fmax=0.05)
