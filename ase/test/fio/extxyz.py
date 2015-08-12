# additional tests of the extended XYZ file I/O
# (which is also included in oi.py test case)
# maintainted by James Kermode <james.kermode@gmail.com>

import os

import numpy as np

import ase.io
from ase.lattice import bulk

# array data of shape (N, 1) squeezed down to shape (N, ) -- bug fixed in commit r4541
at = bulk('Si')
ase.io.write('to.extxyz', at)
at.arrays['ns_extra_data'] = np.zeros( (len(at), 1) )
assert at.arrays["ns_extra_data"].shape == (2, 1)

ase.io.write('to_new.extxyz', at)
at_new = ase.io.read('to_new.extxyz')
assert at_new.arrays["ns_extra_data"].shape == (2, )

os.unlink('to.extxyz')
os.unlink('to_new.extxyz')

# write sequence of images with different numbers of atoms -- bug fixed in commit r4542
images = [at, at*(2, 1, 1), at*(3, 1, 1) ]
ase.io.write('multi.extxyz', images)
read_images = ase.io.read('multi.extxyz@:')
assert read_images == images
os.unlink('multi.extxyz')

# read xyz containing trailing blank line
f = open("structure.extxyz", 'w')
f.write('''4
Coordinates
Mg        -4.25650        3.79180       -2.54123
C         -1.15405        2.86652       -1.26699
C         -5.53758        3.70936        0.63504
C         -7.28250        4.71303       -3.82016

''')
f.close()
a = ase.io.read('structure.extxyz')
os.unlink('structure.extxyz')

# read xyz with / in key value
f = open("slash.extxyz", 'w')
f.write('''4
key1=a key2=a/b
Mg        -4.25650        3.79180       -2.54123
C         -1.15405        2.86652       -1.26699
C         -5.53758        3.70936        0.63504
C         -7.28250        4.71303       -3.82016
''')
f.close()
a = ase.io.read('slash.extxyz')
assert a.info['key1'] == r'a'
assert a.info['key2'] == r'a/b'
os.unlink('slash.extxyz')
