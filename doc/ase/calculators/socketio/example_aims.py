import sys

from ase.build import molecule
from ase.optimize import BFGS
from ase.calculators.aims import Aims
from ase.calculators.socketio import SocketIOCalculator

# Environment-dependent parameters -- please configure according to machine
# Note that FHI-aim support for the i-PI protocol must be specifically
# enabled at compile time, e.g.: make -f Makefile.ipi ipi.mpi
species_dir = '/home/aimsuser/src/fhi-aims.171221_1/species_defaults/light'
command = 'ipi.aims.171221_1.mpi.x'

# This example uses INET; see other examples for how to use UNIX sockets.
port = 31415

atoms = molecule('H2O', vacuum=3.0)
atoms.rattle(stdev=0.1)

aims = Aims(command=command,
            use_pimd_wrapper=('localhost', port),
            compute_forces=True,
            xc='LDA',
            species_dir=species_dir)

opt = BFGS(atoms, trajectory='opt.aims.traj', logfile='opt.aims.log')

with SocketIOCalculator(aims, log=sys.stdout, port=port) as calc:
    atoms.calc = calc
    opt.run(fmax=0.05)
