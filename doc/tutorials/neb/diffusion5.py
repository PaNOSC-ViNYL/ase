from ase.neb import NEBtools
from ase.io import read

images = read('neb.traj@-5:')

nebtools = NEBtools(images)

# Get the calculated barrier and the energy change of the reaction.
Ef, dE = nebtools.get_barrier()

# Get the barrier without any interpolation between highest images.
Ef, dE = nebtools.get_barrier(fit=False)

# Get the actual maximum force at this point in the simulation.
max_force = nebtools.get_fmax()

# Create a figure like that coming from ase-gui.
fig = nebtools.plot_band()
fig.savefig('diffusion-barrier.png')

# Create a figure with custom parameters.
from matplotlib import pyplot, rcParams
rcParams.update({'font.size': 10})
fig = pyplot.figure(figsize=(4.5, 3))
ax = fig.add_axes((0.15, 0.15, 0.8, 0.75))
nebtools.plot_band(ax)
fig.savefig('diffusion-barrier.png')
