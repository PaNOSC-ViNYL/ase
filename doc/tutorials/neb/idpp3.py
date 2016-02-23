import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.calculators.emt import EMT
from ase.neb import NEB
from ase.visualize import view
from ase.optimize.fire import FIRE as QuasiNewton
from ase.lattice.cubic import FaceCenteredCubic

#set the number of images you want
nimages = 5

#some algebra to determine surface normal and the plane of the surface
d3=[2,1,1]
a1=np.array([0,1,1])
d1=np.cross(a1,d3)
a2=np.array([0,-1,1])
d2=np.cross(a2,d3)

#create your slab
slab  =FaceCenteredCubic(directions=[d1,d2,d3],
                         size=(2,1,2),
                         symbol=('Pt'),
                         latticeconstant=3.9)

#add some vacuum to your slab
uc = slab.get_cell()
print(uc)
uc[2] += [0,0,10]  #there are ten layers of vacuum
uc = slab.set_cell(uc,scale_atoms=False)
#view the slab to make sure it is how you expect
view(slab)

#some positions needed to place the atom in the correct place
x1 = 1.379
x2 = 4.137
x3 = 2.759
y1 = 0.0
y2 = 2.238
z1 = 7.165
z2 = 6.439


#Add the adatom to the list of atoms and set constraints of surface atoms.
slab += Atoms('N', [ ((x2+x1)/2,y1,z1+1.5)])
mask = [atom.symbol == 'Pt' for atom in slab]
slab.set_constraint(FixAtoms(mask=mask))

#optimise the initial state
# Atom below step
initial = slab.copy()
initial.set_calculator(EMT())
relax = QuasiNewton(initial)
relax.run(fmax=0.05)
view(initial)

#optimise the initial state
# Atom above step
slab[-1].position = (x3,y2+1,z2+3.5)
final = slab.copy()
final.set_calculator(EMT())
relax = QuasiNewton(final)
relax.run(fmax=0.05)
view(final)

#create a list of images for interpolation
images = [initial]
for i in range(nimages):
    images.append(initial.copy())

for image in images:
    image.set_calculator(EMT())

images.append(final)
view(images)

#carry out idpp interpolation
neb = NEB(images)
neb.interpolate('idpp')

#Run NEB calculation
qn = QuasiNewton(neb, trajectory='N_diffusion.traj',  logfile='N_diffusion.log')
qn.run(fmax=0.05)


