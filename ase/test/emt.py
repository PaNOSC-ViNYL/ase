import numpy as npy
from ase import *

a = 3.60
b = a / 2
cu = Atoms(positions=[(0, 0, 0)],
           symbols='Cu',
           cell=[(0, b, b),
                 (b, 0, b),
                 (b, b, 0)],
           pbc=1,
           calculator=EMT())
e0 = cu.get_potential_energy()
print e0

cu.set_cell(cu.get_cell() * 1.001)
e1 = cu.get_potential_energy()
V = a**3 / 4
B = 2 * (e1 - e0) / 0.003**2 / V * 160.2
print B

for i in range(4):
    x = 0.001 * i
    A = npy.array([(x, b, b+x),
                   (b, 0, b),
                   (b, b, 0)])
    cu.set_cell(A)
    e = cu.get_potential_energy() - e0
    if i == 0:
        print i, e
    else:
        print i, e, e / x**2

A = npy.array([(0, b, b),
               (b, 0, b),
               (6*b, 6*b, 0)])
R = npy.zeros((2, 3))
for i in range(1, 2):
    R[i] = i * A[2] / 6
print (Atoms(positions=R, symbols='Cu2',
             pbc=1, cell=A, calculator=EMT()).get_potential_energy() - 2 * e0) / 2

A = npy.array([(0, b, b),
               (b, 0, b),
               (10*b, 10*b, 0)])
R = npy.zeros((3, 3))
for i in range(1, 3):
    R[i] = i * A[2] / 10
print (Atoms(positions=R, symbols='Cu3',
             pbc=1, cell=A, calculator=EMT()).get_potential_energy() - 3 * e0) / 2

A = npy.array([(0, b, b),
               (b, 0, b),
               (b, b, 0)])
R = npy.zeros((3, 3))
for i in range(1, 3):
    R[i] = i * A[2]
print (Atoms(positions=R, symbols='Cu3',
             pbc=(1, 1, 0), cell=A, calculator=EMT()).get_potential_energy() - 3 * e0) / 2

