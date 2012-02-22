import os

from ase.test import NotAvailable

try:
    nwchem_basis_library = os.getenv('NWCHEM_BASIS_LIBRARY')
    if nwchem_basis_library == None:
        raise NotAvailable('NWCHEM_BASIS_LIBRARY not defined')
    nwchem_command = os.getenv('NWCHEM_COMMAND')
    if nwchem_command == None:
        raise NotAvailable('NWCHEM_COMMAND not defined')
except NotAvailable:
    raise NotAvailable('Nwchem required')

import numpy as np

from ase.tasks.main import run

atoms, task = run("nwchem molecule O2 O -l -p task=gradient")
atoms, task = run('nwchem molecule O2 O -s')
ae = 2 * np.array(task.results['O'])[0] - np.array(task.results['O2'])[0]
assert abs(ae - 6.605) < 1e-3
