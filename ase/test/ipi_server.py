import sys
import threading
from subprocess import Popen

import numpy as np

from ase.calculators.ipi import IPIClient, IPICalculator
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from ase.cluster.icosahedron import Icosahedron


def getatoms():
    return Icosahedron('Au', 3)


def run_server(launchclient=True):
    atoms = getatoms()

    with IPICalculator(log=sys.stdout) as calc:
        if launchclient:
            thread = launch_client_thread()
        atoms.calc = calc
        opt = BFGS(atoms)
        opt.run()

    if launchclient:
        thread.join()

    forces = atoms.get_forces()
    energy = atoms.get_potential_energy()

    atoms.calc = EMT()
    ref_forces = atoms.get_forces()
    ref_energy = atoms.get_potential_energy()

    refatoms = run_normal()
    ref_energy = refatoms.get_potential_energy()
    eerr = abs(energy - ref_energy)
    ferr = np.abs(forces - ref_forces).max()

    perr = np.abs(refatoms.positions - atoms.positions).max()
    print('errs e={} f={} pos={}'.format(eerr, ferr, perr))
    assert eerr < 1e-12, eerr
    assert ferr < 1e-12, ferr
    assert perr < 1e-12, perr

def run_normal():
    atoms = getatoms()
    atoms.calc = EMT()
    opt = BFGS(atoms)
    opt.run()
    return atoms

def run_client():
    import time
    time.sleep(1)
    atoms = getatoms()
    atoms.calc = EMT()
    with open('client.log', 'w') as fd:
        client = IPIClient(log=fd)
        client.run(atoms, use_stress=False)


def launch_client_thread():
    thread = threading.Thread(target=run_client)
    thread.start()
    return thread


import sys
if len(sys.argv) > 1:
    arg = sys.argv[1]
    if arg == 'server':
        run_server(False)
    elif arg == 'client':
        run_client()
    else:
        assert 0
else:
    run_server()
