from ase.build import molecule
from ase.neb import NEB
from ase.calculators.emt import EMT
from ase.optimize.fire import FIRE as QuasiNewton
from ase.visualize import view

#Optimise molecule
initial = molecule('C2H6')
initial.set_calculator(EMT())
relax = QuasiNewton(initial)
relax.run(fmax=0.05)
view(initial)

#Create final state
final = initial.copy()
final.positions[2:5] = initial.positions[[3, 4, 2]]

#Generate blank images
images = [initial]

for i in range(9):
    images.append(initial.copy())

for image in images:
    image.set_calculator(EMT())
   
images.append(final)

#Run IDPP interpolation
neb = NEB(images)
neb.interpolate('idpp')

#Run NEB calculation
qn = QuasiNewton(neb, trajectory='ethane_idpp.traj',  logfile='ethane_idpp.log')
qn.run(fmax=0.05)
