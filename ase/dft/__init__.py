import numpy as npy
from ase.dft.stm import STM
from ase.dft.dos import DOS
from ase.dft.wannier import Wannier


def monkhorst_pack(size):
    if npy.less_equal(size, 0).any():
        raise ValueError('Illegal size: %s' % list(size))
    kpts = npy.indices(size).transpose((1, 2, 3, 0)).reshape((-1, 3))
    return (kpts + 0.5) / size - 0.5

def get_distribution_moment(x, y, order=0):
    """Return the moment of nth order of distribution.
    
    1st and 2nd order moments of a band correspond to the band's
    center and width respectively.
    
    For integration, the trapezoid rule is used.
    """

    x = npy.array(x)
    y = npy.array(y)
        
    if order==0:
        return npy.trapz(y, x)
    elif type(order) is int:
        return npy.trapz(x**order * y, x) / npy.trapz(y, x)
    elif hasattr(order, '__iter__'):
        return [get_distribution_moment(x, y, n) for n in order]
    else:
        raise ValueError('Illegal order: %s' % str(order))

