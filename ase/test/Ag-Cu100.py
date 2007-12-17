from ase import *

# Distance between Cu atoms on a (100) surface:
d = 3.6 / sqrt(2)
initial = Atoms(symbols='Cu',
                positions=[(0, 0, 0)],
                cell=(d, d, 10),
                pbc=(True, True, False))
initial *= (5, 5, 1)  # 5x5 (100) surface-cell

# Approximate height of Ag atom on Cu(100) surfece:
h0 = 2.0
initial += Atom('Ag', (d / 2, d / 2, h0))

if 0:
    view(initial)

# Make band:
images = [initial.copy() for i in range(6)]
neb = NEB(images, climb=True)

# Set constraints and calculator:
constraint = FixAtoms(range(len(initial) - 1))
for image in images:
    image.set_calculator(ASAP())
    image.constraints.append(constraint)

# Displace last image:
images[-1].positions[-1] += (d, 0, 0)

# Relax height of Ag atom for initial and final states:
for image in [images[0], images[-1]]:
    QuasiNewton(image).run(fmax=0.01)

# Interpolate positions between initial and final states:
neb.interpolate()

for image in images:
    print image.positions[-1], image.get_potential_energy()

traj = PickleTrajectory('mep.traj', 'w')

dyn = MDMin(neb, dt=0.1)
dyn = QuasiNewton(neb)
#dyn.attach(neb.writer(traj))
dyn.run(fmax=0.01, steps=40)

for image in images:
    print image.positions[-1], image.get_potential_energy()
