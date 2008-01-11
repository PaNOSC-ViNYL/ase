import os
import pickle
import tempfile

import numpy as npy

from ase.io.trajectory import write_trajectory


def gui(atoms):
    fd, filename = tempfile.mkstemp('.traj', 'ag-')
    os.close(fd)
    if not isinstance(atoms, list):
        atoms = [atoms]
    write_trajectory(filename, atoms)
    os.system('(ag %s &); (sleep 5; rm %s) &' % (filename, filename))
