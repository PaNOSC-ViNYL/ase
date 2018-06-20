import sys
from subprocess import Popen

from ase.build import molecule
from ase.io import write
from ase.optimize import BFGS
from ase.calculators.socketio import SocketIOCalculator

unixsocket = 'ase_server_socket'
targetscript = 'example_client_gpaw.py'

atoms = molecule('H2O', vacuum=3.0)
atoms.rattle(stdev=0.1)
write('initial.traj', atoms)

opt = BFGS(atoms, trajectory='opt.driver.traj', logfile='opt.driver.log')

with SocketIOCalculator(log=sys.stdout,
                        unixsocket=unixsocket) as calc:
    # Server is now running and waiting for connections.  Now is a
    # good time to launch the client process.  Uncomment the
    # following line to do so:
    # proc = Popen([sys.executable, targetscript])
    #
    # Else the server will wait for the user to run
    # that file (or another client) manually.

    atoms.calc = calc
    opt.run(fmax=0.05)
