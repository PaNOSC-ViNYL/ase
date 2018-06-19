import sys
from ase.build import molecule
from ase.calculators.siesta import Siesta
from ase.optimize import BFGS
from ase.calculators.ipi import IPICalculator

unixsocket = 'siesta'

fdf_arguments = {'MD.TypeOfRun': 'Master',
                 'Master.code': 'i-pi',
                 'Master.interface': 'socket',
                 'Master.address': unixsocket,
                 #'Master.port': port,
                 'Master.socketType': 'unix'}

atoms = molecule('H2O', vacuum=3.0)
atoms.rattle(stdev=0.1)

siesta = Siesta(fdf_arguments=fdf_arguments)
opt = BFGS(atoms, trajectory='opt.siesta.traj', logfile='opt.siesta.log')

with IPICalculator(siesta, log=sys.stdout, unixsocket=unixsocket) as calc:
    atoms.calc = calc
    opt.run(fmax=0.05)
