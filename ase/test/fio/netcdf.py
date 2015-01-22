import numpy as np

from ase.test import NotAvailable

try:
    from scipy.io.netcdf import NetCDFFile
except ImportError:
    raise NotAvailable('Scipy too old')

# Write array
a1 = np.random.rand(5, 5)
a2 = a1 * 2 - 5
nc = NetCDFFile('test.nc', 'w')
nc.createDimension('dimx', a1.shape[0])
nc.createDimension('dimy', a1.shape[1])
nc.createVariable('matrix1', 'd', ('dimx', 'dimy'))[:] = a1
nc.createVariable('matrix2', 'd', ('dimx', 'dimy'))[:] = a2
nc.sync()
nc.close()

# Read array
nc = NetCDFFile('test.nc', 'r')
b1 = nc.variables['matrix1'][:]
b2 = nc.variables['matrix2'][:]

assert np.all(a1 == b1) and np.all(a2 == b2)
