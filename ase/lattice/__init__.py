import functools
from ase.utils.deprecate import deprecate
from ase.build.bulk import bulk as newbulk
__all__ = ['bulk']


@functools.wraps(newbulk)
def bulk(*args, **kwargs):
    deprecate('Use ase.build.bulk() instead', '3.11')
    return newbulk(*args, **kwargs)
