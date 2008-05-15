.. module:: io
   :synopsis: File input-output module


File input and output
=====================

The ``io`` module has two basic methods: ``read`` and ``write``, both
of which are accessed by either::

  >>> from ase import *

or::
  
  >>> from ase.io import read, write

The two methods are described below:

.. function:: read(filename, index=-1, format=None)
    
  Read Atoms object(s) from file.

  ::

    filename: str
        Name of the file to read from.
    index: int or slice
        If the file contains several configurations, the last configuration
        will be returned by default.  Use index=n to get configuration
        number n (counting from zero).
    format: str
        Used to specify the file-format.  If not given, the file-format
        will be guessed.

  The accepted input formats:

  =========================  ===========
  format                     short name
  =========================  ===========
  GPAW restart-file          gpw
  Dacapo netCDF output file  dacapo
  Old ASE netCDF trajectory  nc
  Virtual Nano Lab file      vnl
  ASE pickle trajectory      traj
  GPAW text output           gpaw-text
  CUBE file                  cube
  Dacapo text output         dacapo-text
  XYZ-file                   xyz
  =========================  ===========


The ``read`` method is only designed to retrive the atom configuration
from a file, but for the `cube` format you can import read_cube::

  >>> from ase.io.cube import read_cube

which takes the additional keyword flag ``read_data`` which if set to
``True`` causes ``read_cube`` to return a ``(data, atoms)`` tuple.

.. function:: write(filename, images, format=None, **kwargs)
   
  Write Atoms object(s) to file.

  ::

    filename: str
        Name of the file to write to.
    images: Atoms object or list of Atoms objects
        A single Atoms object or a list of Atoms objects.
    format: str
        Used to specify the file-format.  If not given, the file-format
        will be taken from suffix of the filename.

  The accepted output formats:

  =========================  ===========
  format                     short name
  =========================  ===========
  Old ASE netCDF trajectory  nc
  ASE pickle trajectory      traj
  CUBE file                  cube
  XYZ-file                   xyz
  Protein Data Bank          pdb
  gOpenMol .plt file         plt  
  Python script              py

  Encapsulated Postscript    eps
  Portable Network Graphics  png
  Persistance of Vision      pov
  =========================  ===========

  The use of additional keywords is format specific.

  The ``cube`` and ``plt`` formats accept (plt requires it) a ``data``
  keyword, which can be used to write a 3D array to the file along
  with the nuclei coordinates. The array must be real-valued.

  The ``eps``, ``png``, and ``pov`` formats are all graphics formats,
  and accept the additional keywords::

    rotation='', show_unit_cell=0, radii=None, bbox=None, colors=None

  ::

    rotation: str
      The rotation angles, e.g. '45x,70y,90z'
    show_unit_cell: int
      Can be 0, 1, 2 to either not show, show, or show all of the unit cell
    radii: array
      An array of same length as the list of atoms, indicating the sphere radii
    bbox: array
      XXX
    colors: array
      An array of same length as the list of atoms, indicating the rgb color
      code for each atom

  The ``pov`` accepts the additional keywords:
  
  XXX

  For ``pov`` the elements of the color array can also be strings, or 4,
  and 5 vectors.

  XXX

