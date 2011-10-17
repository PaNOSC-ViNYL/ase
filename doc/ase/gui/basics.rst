.. module:: basics

==================================
ag basics and command line options
==================================

General use
-----------

Visualizing a system with ag is straight-forward using a regular
mouse. The scroll function allows to change the magnification, the
left mouse button selects atoms, the right mouse button allows to
rotate, and the middle button allows to translate the system on the
screen. 

Depending on the number of selected atoms, ag automatically measures
different quantities: 

================================= ======================================
Selection			  measurement
================================= ======================================
single atom                       xyz position and atomic symbol
two atoms                         interatomic distance and symbols
three atoms                       all three internal angles and
      				  symbols 
four atoms, selected sequentially Measures the dihedral angle,
     	    	     		  e.g. the angle between bonds 12 and 34
more than four atoms		  chemical composition of selection. 
================================= ======================================

ag can save the following file formats: 

=========== =================================
File format Comment
=========== =================================
xyz 	    XYZ file
traj	    ASE trajectory
pdb	    PDB file
cube	    Gaussian cube file
py 	    Python script
vnl	    VNL file
png	    Portable Network Graphics
pov	    Persistance of Vision
eps	    Encapsulated PostScript
in	    FHI-aims geometry input
POSCAR	    VASP geometry input
bundle	    ASE bundle trajectory
cif	    Crystallographic Information File
=========== =================================

Files
-----

The :program:`ag` program can read all the file formats the ASE's
:func:`~ase.io.read` function can understand.

::
  
  $ ag N2Fe110-path.traj


Selecting part of a trajectory
------------------------------
  
A Python-like syntax for selecting a subset of configurations can be
used.  Instead of the Python syntax ``list[start:stop:step]``, you use
:file:`filaname@start:stop:step`::

  $ ag x.traj@0:10:1  # first 10 images
  $ ag x.traj@0:10    # first 10 images
  $ ag x.traj@:10     # first 10 images
  $ ag x.traj@-10:    # last 10 images
  $ ag x.traj@0       # first image
  $ ag x.traj@-1      # last image
  $ ag x.traj@::2     # every second image

If you want to select the same range from many files, the you can use
the :option:`-n` or :option:`--image-number` option::

  $ ag -n -1 *.traj   # last image from all files
  $ ag -n 0 *.traj    # first image from all files

.. tip::

  Type :program:`ag -h` for a description of all command line options.


Writing files
-------------

::

  $ ag -n -1 a*.traj -o new.traj

Possible formats are: ``traj``, ``xyz``, ``cube``, ``pdb``, ``eps``,
``png``, and ``pov``.  For details, see the :mod:`~ase.io` module
documentation.

Interactive use
---------------

The :program:`ag` program can also be launched directly from a Python
script or interactive session:

>>> from ase import *
>>> atoms = ...
>>> view(atoms)

or

>>> view(atoms, repeat=(3, 3, 2))

or, to keep changes to your atoms:

>>> atoms.edit()


NEB calculations
----------------

Use :menuselection:`Tools --> NEB` to plot energy barrier.

::
  
  $ ag --interpolate 3 initial.xyz final.xyz -o interpolated_path.traj


Plotting data from the command line
-----------------------------------
Plot the energy relative to the energy of the first image as a
function of the distance between atom 0 and 5::

  $ ag -g "d(0,5),e-E[0]" x.traj
  $ ag -t -g "d(0,5),e-E[0]" x.traj > x.dat  # No GUI, write data to stdout

The symbols are the same as used in the plotting data function. 


Defaults for ag
---------------

Using a file ~/.ase/gui.py, certain defaults can be set. If it exists,
this file is executed after initializing the variables and colours
normally used in ag. One can change the default graphs that are
plotted, and the default radii for displaying specific atoms. This
example will display the energy evolution and the maximal force in a
graph and also display Cu atoms (Z=29) with a radius of 1.6 Angstrom.

::

  gui_default_settings['gui_graphs_string'] = "i, e - min(E), fmax"
  gui_default_settings['covalent_radii'] = [[29,1.6]]
