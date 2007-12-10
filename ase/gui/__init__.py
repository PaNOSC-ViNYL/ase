import os
import pickle
import tempfile

import numpy as npy

def gui(atoms):
    fd, filename = tempfile.mkstemp('.pckl', 'ag-')
    if not isinstance(atoms, list):
        atoms = [atoms]
    os.write(fd, pickle.dumps(atoms))
    os.close(fd)
    os.system('(ag --read-pickled-data-from-file %s &); (sleep 5; rm %s) &' %
              (filename, filename))
