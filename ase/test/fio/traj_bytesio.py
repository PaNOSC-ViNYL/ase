import io
from ase.build import bulk
from ase.collections import g2
from ase.io import iread, write
from ase.io.trajectory import bytestoimages, imagestobytes

images = [bulk('Si') + bulk('Fe')] + list(g2)

txt = imagestobytes(images)
images2 = bytestoimages(txt)

for atoms1, atoms2 in zip(images, images2):
    assert atoms1 == atoms2
