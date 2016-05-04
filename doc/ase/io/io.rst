.. module:: ase.io
   :synopsis: File input-output module

=====================
File input and output
=====================

.. seealso::
    
    * :mod:`ase.io.trajectory`
    
.. toctree::
    :hidden:
    
    trajectory
    opls
    
    
The :mod:`ase.io` module has three basic functions: :func:`read`,
:func:`iread` and :func:`write`. The methods are described here:

.. autofunction:: read
.. autofunction:: iread
.. autofunction:: write

These are the file-formats that are recognized (formats with a ``+`` support
multiple configurations):

.. csv-table::
    :file: io.csv
    :header-rows: 1
    
.. note::

    Even though that ASE does a good job reading the above listed
    formats, it may not read some unusual features or strangely
    formatted files.

    For the CIF format, STAR extensions as save frames, global blocks,
    nested loops and multi-data values are not supported.  Furthermore,
    ASE currently assumes the ``loop_`` identifier, and the following
    loop variable names to be on separate lines.

The :func:`read` function is only designed to retrieve the atomic configuration
from a file, but for the CUBE format you can import the function:

.. function:: read_cube_data


which will return a ``(data, atoms)`` tuple::

  from ase.io.cube import read_cube_data
  data, atoms = read_cube_data('abc.cube')


Examples
========

::

    from ase.build import *
    adsorbate = Atoms('CO')
    adsorbate[1].z = 1.1
    a = 3.61
    slab = fcc111('Cu', (2, 2, 3), a=a, vacuum=7.0)
    add_adsorbate(slab, adsorbate, 1.8, 'ontop')

Write PNG image::

    write('slab.png', slab * (3, 3, 1), rotation='10z,-80x')

.. image:: io1.png

Write POVRAY file::

    write('slab.pov', slab * (3, 3, 1), rotation='10z,-80x')

This will write both a ``slab.pov`` and a ``slab.ini`` file.  Convert
to PNG with the command ``povray slab.ini`` or use the
``run_povray=True`` option:

.. image:: io2.png

Here is an example using ``bbox``::

    d = a / 2**0.5
    write('slab.pov', slab * (2, 2, 1),
          bbox=(d, 0, 3 * d, d * 3**0.5))

.. image:: io3.png

Note that the XYZ-format does not contain information about the unic cell:

>>> write('slab.xyz', slab)
>>> a = read('slab.xyz')
>>> a.get_cell()
array([[ 1.,  0.,  0.],
       [ 0.,  1.,  0.],
       [ 0.,  0.,  1.]])
>>> a.get_pbc()
array([False, False, False], dtype=bool)

Use ASE's native format for writing all information:

>>> write('slab.traj', slab)
>>> b = read('slab.traj')
>>> b.cell
array([[  5.10531096e+00,  -4.11836034e-16,   1.99569088e-16],
       [  2.55265548e+00,   4.42132899e+00,   7.11236625e-17],
       [  8.11559027e+00,   4.68553823e+00,   1.32527034e+01]])
>>> b.pbc
array([ True,  True,  True], dtype=bool)

A script showing all of the povray parameters, and generating the image below,
can be found here: :download:`save_pov.py`

.. image:: NaCl_C6H6.png

An other example showing how to change colors and textures in pov can
be found here: :download:`../../tutorials/saving_graphics.py`.


Adding a new file-format to ASE
===============================

Try to model the read/write functions after the *xyz* format as implemented
in :git:`ase/io/xyz.py` and also read, understand and update
:git:`ase/io/formats.py`.
