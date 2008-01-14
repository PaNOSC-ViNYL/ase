import numpy as npy
from ase.dft.stm import STM
from ase.dft.dos import DOS


def monkhorst_pack(size):
    if npy.less_equal(size, 0).any():
        raise ValueError('Illegal size: %s' % list(size))
    kpts = npy.indices(size).transpose((1, 2, 3, 0)).reshape((-1, 3))
    return (kpts + 0.5) / size - 0.5
