def write_plt(filename, array, cell=None):
    cell = npy.asarray(cell, float)
    if cell.ndim == 2:
        c = cell.copy()
        cell = c.diagonal()
        c.flat[::4] = 0.0
        assert not c.any()
    dims = npy.array(array.shape)
    f = open(filename, 'w')
    npy.array([3, 4]).tofile(f)
    dims.tofile(f)
