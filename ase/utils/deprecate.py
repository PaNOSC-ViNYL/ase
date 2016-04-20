import warnings
from distutils.version import StrictVersion

from ase import __version__


def deprecate(msg, version):
    if 1:#StrictVersion(__version__) >= version:
        warnings.warn(msg, stacklevel=2)
