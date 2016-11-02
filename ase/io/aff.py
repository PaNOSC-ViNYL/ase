import warnings
from ase.io.ulm import (ulmopen as affopen,
                        InvalidULMFileError as InvalidAFFError,
                        Reader, Writer, DummyWriter)

__all__ = ['affopen', 'InvalidAFFError',
           'Reader', 'Writer', 'DummyWriter']

warnings.warn('Please use ase.io.ulm instead.')
