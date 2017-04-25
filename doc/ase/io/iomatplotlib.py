# creates: iomatplotlib1.png, iomatplotlib2.png

import matplotlib.pyplot as plt
from ase.io.matplotlib import write_matplotlib
from ase.lattice.cubic import FaceCenteredCubic

slab = FaceCenteredCubic('Au', size=(2, 2, 2))
fig, ax = plt.subplots()
write_matplotlib(slab, ax, radii=0.3, rotation=('90x,45y,0z'))
fig.savefig("iomatplotlib1.png")

slab = FaceCenteredCubic('Au', size=(2, 2, 2))
fig, axarr = plt.subplots(1, 4, figsize=(15, 5))
write_matplotlib(slab, axarr[0], radii=0.3, rotation=('0x,0y,0z'))
write_matplotlib(slab, axarr[1], radii=0.3, rotation=('45x,0y,0z'))
write_matplotlib(slab, axarr[2], radii=0.3, rotation=('45x,45y,0z'))
write_matplotlib(slab, axarr[3], radii=0.3, rotation=('0x,0y,0z'))
axarr[0].set_title("No rotation")
axarr[1].set_xlabel("X-axis, [$\mathrm{\AA}$]")
axarr[1].set_ylabel("Y-axis, [$\mathrm{\AA}$]")
axarr[2].set_axis_off()
axarr[3].set_xlim(2, 6)
axarr[3].set_ylim(2, 6)
fig.savefig("iomatplotlib2.png")
