from ase.dft import get_distribution_moment
import numpy as npy

precision = 1E-8

x = npy.linspace(-50., 50., 1000)
y = map(lambda x: npy.exp(-x**2 / 2.), x)
area, center, mom2 = get_distribution_moment(x, y, (0, 1, 2))
assert sum((abs(area - npy.sqrt(2. * npy.pi)), abs(center), abs(mom2 - 1.))) < precision

x = npy.linspace(-1., 1., 100000)
for order in range(0, 9):
    y = map(lambda x: x**order, x)
    area = get_distribution_moment(x, y)
    assert abs(area - (1. - (-1.)**(order + 1)) / (order + 1.)) < precision

x = npy.linspace(-50., 50., 100)
y = map(lambda x: npy.exp(-2. * (x - 7.)**2 / 10.) + npy.exp(-2. * (x + 5.)**2 / 10.), x)
center=get_distribution_moment(x, y, 1)
assert abs(center - 1.) < precision 
