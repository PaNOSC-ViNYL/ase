import sys
from ase.build import molecule
from ase.calculators.aims import Aims
from ase.optimize import BFGS
from ase.calculators.ipi import IPICalculator

port = 27182
species_dir = '/home/aimsuser/src/fhi-aims.171221_1/species_defaults/light'
command = 'ipi.aims.171221_1.mpi.x',

atoms = molecule('H2O', vacuum=3.0)
atoms.rattle(stdev=0.1)

calc = Aims(command=command,
            use_pimd_wrapper=('localhost', port),
            compute_forces=True,
            xc='LDA',
            species_dir=species_dir)

opt = BFGS(atoms, trajectory='opt.aims.traj', logfile='opt.aims.log')

with IPICalculator(calc, log=sys.stdout, port=port) as calc:
    atoms.calc = calc
    opt.run(fmax=0.05)
