import os
import pickle
import tempfile

import numpy as npy

from ase.io.trajectory import write_trajectory


def gui(atoms):
    if not isinstance(atoms, list):
        atoms = [atoms]
    filename = tempfile.mktemp('.traj', 'ag-')
    write_trajectory(filename, atoms)
    os.system('(ag %s &); (sleep 15; rm %s) &' % (filename, filename))
