import os
import pickle
import tempfile

import numpy as npy

from ase.io.trajectory import write_trajectory


def gui(atoms):
    if not isinstance(atoms, list):
        atoms = [atoms]
    filename = tempfile.mktemp('.traj', 'ag-')
    calc = atoms.get_calculator()
    atoms.set_calculator(None)
    write_trajectory(filename, atoms)
    atoms.set_calculator(calc)
    os.system('(ag %s &); (sleep 15; rm %s) &' % (filename, filename))
