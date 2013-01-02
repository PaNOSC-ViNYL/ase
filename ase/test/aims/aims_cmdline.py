from ase.test import NotAvailable

from ase.calculators.aims import Aims

if Aims().get_command() is None:
    raise NotAvailable('FHI-aims required')

from ase.tasks.main import run

from os import environ

# warning! parameters are not converged - only an illustration!

atoms, task = run('aims bulk Li -x bcc -a 3.6 --k-point-density 1.5 --srelax 0.05 --srelaxsteps 1 -t stress -p xc=pw-lda,sc_accuracy_eev=1.e-3,relativistic=none,compute_analytical_stress=True')
atoms, task = run('aims bulk Li -x bcc -a 3.6 -t stress -s')
data = task.data['Li']
