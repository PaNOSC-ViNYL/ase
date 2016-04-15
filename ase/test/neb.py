import threading

from ase.test import World
from ase.io import Trajectory, read
from ase.neb import NEB, NEBtools
from ase.calculators.morse import MorsePotential
from ase.calculators.singlepoint import SinglePointCalculator
from ase.optimize import BFGS


fmax = 0.05
nimages = 3

print([a.get_potential_energy() for a in Trajectory('H.traj')])
images = [Trajectory('H.traj')[-1]]
for i in range(nimages):
    images.append(images[0].copy())
images[-1].positions[6, 1] = 2 - images[0].positions[6, 1]
neb = NEB(images)
neb.interpolate()
if 0:  # verify that initial images make sense
    from ase.visualize import view
    view(neb.images)

for image in images:
    image.set_calculator(MorsePotential())

dyn = BFGS(neb, trajectory='mep.traj')  # , logfile='mep.log')

dyn.run(fmax=fmax)

for a in neb.images:
    print(a.positions[-1], a.get_potential_energy())

results = [neb.emax]

neb.climb = True
dyn.run(fmax=fmax)
results.append(neb.emax)

# Check NEB tools.
nt_images = read('mep.traj@-4:')
nebtools = NEBtools(nt_images)
nt_fmax = nebtools.get_fmax(climb=True)
Ef, dE = nebtools.get_barrier()
print(Ef, dE, fmax, nt_fmax)
assert nt_fmax < fmax
assert abs(Ef - 1.389) < 0.001


# Test NEB in parallel using some tricks and two threads:
def run_neb_calculation(cpu):
    images = [Trajectory('H.traj')[-1]]
    for i in range(nimages):
        images.append(images[0].copy())
    images[-1].positions[6, 1] = 2 - images[0].positions[6, 1]
    images[-1].set_calculator(
        SinglePointCalculator(energy=images[0].get_potential_energy(),
                              atoms=images[-1]))
    neb = NEB(images, parallel=True, world=cpu)
    neb.interpolate()

    images[cpu.rank + 1].set_calculator(MorsePotential())

    dyn = BFGS(neb)
    dyn.run(fmax=fmax)

    if cpu.rank == 1:
        results.append(neb.emax)
    
    neb.climb = True
    dyn.run(fmax=fmax)
    if cpu.rank == 1:
        results.append(neb.emax)
    
w = World(nimages - 1)
ranks = [w.get_rank(r) for r in range(w.size)]
threads = [threading.Thread(target=run_neb_calculation, args=(rank,))
           for rank in ranks]
for t in threads:
    t.start()
for t in threads:
    t.join()

print(results)
assert abs(results[0] - results[2]) < 1e-5
assert abs(results[1] - results[3]) < 1e-5
