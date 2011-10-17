import numpy as np
from ase.tasks.main import run
atoms, task = run('H2 -R 0.01 -F --atomize')
atoms, task = run('H2 H -s')
results = np.array(task.results['H2'])
assert abs(results -
           [1.1589, 0.0882, 0.7789, 869.0545, 5.2611, 5.3493]).max() < 0.001

atoms, task = run('bulk Cu -F')
atoms, task = run('bulk Cu -s')
results = np.array(task.results['Cu'])
assert abs(results -
           [-0.0057, 0.0014, 11.5654, 134.4389]).max() < 0.001
