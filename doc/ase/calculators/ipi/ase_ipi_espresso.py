from ase.build import molecule
from ase.optimize import BFGS
from ase.calculators.espresso import Espresso

atoms = molecule('H2O', vacuum=3.0)
atoms.rattle(stdev=0.1)

# Environment-dependent parameters:
pseudopotentials = {'H': 'H.pbe-rrkjus.UPF',
                    'O': 'O.pbe-rrkjus.UPF'}
command = 'pw.x < {prefix}.pwi -ipi {host}:{port} > {prefix}.pwo'
pseudo_dir='.'

espresso = Espresso(command=command,
                    ecutwfc=30.0,
                    pseudopotentials=pseudopotentials,
                    pseudo_dir=pseudo_dir)

opt = BFGS(atoms, trajectory='opt.traj',
           logfile='opt.log')

with espresso.ipi() as calc:
    atoms.calc = calc
    opt.run(fmax=0.05)

# Note: As of now, the QE process exits with status 128
# even if apparently successful.
