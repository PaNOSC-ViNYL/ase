# creates: Al_phonon.png, Al_mode.gif
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.dft.kpoints import bandpath
from ase.phonons import Phonons

# Setup crystal and EMT calculator
atoms = bulk('Al', 'fcc', a=4.05)

# Phonon calculator
N = 7
ph = Phonons(atoms, EMT(), supercell=(N, N, N), delta=0.05)
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True)
ph.clean()

path = bandpath('GXULGK', atoms.cell, 100)[0]
bs = ph.get_band_structure(path)

dos = ph.get_dos(kpts=(20, 20, 20)).sample_grid(npts=100, width=1e-3)

# Plot the band structure and DOS:
import matplotlib.pyplot as plt
fig = plt.figure(1, figsize=(7, 4))
ax = fig.add_axes([.12, .07, .67, .85])

emax = 0.035
bs.plot(ax=ax, show=False, emin=0.0, emax=emax)

dosax = fig.add_axes([.8, .07, .17, .85])
dosax.fill_between(dos.weights[0], dos.energy, y2=0, color='grey',
                   edgecolor='k', lw=1)

dosax.set_ylim(0, emax)
dosax.set_yticks([])
dosax.set_xticks([])
dosax.set_xlabel("DOS", fontsize=18)

fig.savefig('Al_phonon.png')

# --- End of literalinclude ---

# Literalinclude start modes
# Write modes for specific q-vector to trajectory files
from ase.dft.kpoints import ibz_points
L = ibz_points['fcc']['L']
ph.write_modes([l / 2 for l in L], branches=[2], repeat=(8, 8, 8), kT=3e-4,
               center=True)
# Literalinclude end modes


# Generate png animation (requires the convert command from ImageMagick):
from subprocess import call
from ase.io.trajectory import Trajectory
from ase.io import write
trajectory = Trajectory('phonon.mode.2.traj', 'r')

for i, atoms in enumerate(trajectory):
    write('picture%02i.png' % i, atoms, show_unit_cell=2,
          rotation='-36x,26.5y,-25z')
    # Flatten images for better quality
    call(['convert', '-flatten', 'picture%02i.png' % i, 'picture%02i.png' % i])

# Concatenate to gif animation
call(['convert', '-delay', '5', '-loop', '0', '-dispose', 'Previous',
      'picture*.png', 'Al_mode.gif'])
