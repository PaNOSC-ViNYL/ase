import pathlib as pl
import sys
from ase.utils import devnull
from ase.parallel import paropen


def is_file_like(fd):
    """Test is a function is already a file-type object.
    Function borrowed from Pandas"""
    if not (hasattr(fd, 'read') or hasattr(fd, 'write')):
        return False

    if not hasattr(fd, '__iter__'):
        return False

    return True


def pathify(fd):
    """Convert none file-like objects to Path objects.
    If fd is a like-like object, it is just returned."""

    if is_file_like(fd):
        # fd is already a buffer
        return fd

    if isinstance(fd, pl.PurePath):
        # fd is already a pathlib object
        return fd

    return pl.Path(fd)


def open_fd(name, mode='w'):
    if name == '-':
        name = sys.stdout
    else:
        name = pathify(name)
    return paropen(name, mode=mode)


class ASEbuffer:
    def __init__ (self, filename, mode='w'):
    """Class for keeping track of internally opened files.
    Should be used as a context manager.
    Only closes filename if we opened file internally."""
        self.mode = mode
        self.opened_internally = False
        self.filename = pathify(filename)

        self.open()

    def open(self):
        if is_file_like(self.filename):
            buffer = self.filename
        else:
            self._opened_internally = True
            buffer = self.filename.open(self.mode)
        self.buffer = buffer
        self.closed = self.buffer.closed
        return self.buffer

    def close(self):
        self.buffer.close()
        self.closed = self.buffer.closed

    def write(self, text):
        self.buffer.write(text)

    def __enter__ (self):
        return self.open()

    def __exit__ (self, exc_type, exc_value, traceback):
        """Only close buffer is we opened it"""
        if self.opened_internally:
            self.close()
