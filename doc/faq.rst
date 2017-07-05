.. _faq:

==========================
Frequently Asked Questions
==========================


ASE-GUI
=======

See also the :mod:`ase.gui`.


How do I export images from a trajectory to png or pov files?
-------------------------------------------------------------

With ase-gui, you can choose :menuselection:`File --> Save`, but this is
not fun if you need to do it for many images.  Here is how to do it on
the command line for a number of images::

  ase gui images.traj@0 -o image0.pov
  ase gui images.traj@1 -o image1.pov
  ase gui images.traj@2 -o image2.pov

If you have many images, it will be easier to do it using the Python
interpreter:

>>> from ase.io import read, write
>>> for n, image in enumerate(read('images.traj@:3')):
...     write('image%d.pov' % n, image, run_povray=True, pause=False,
...           rotation='-90x,10z')

Here, we also:

* run povray to generate png files

* disable pausing between the images

* set a rotation (choose :menuselection:`View --> Rotate ...` in ase-gui to
  select the best rotation angles)

Try:

>>> help(write)

to see all possibilities or read more :func:`here <ase.io.write>`.



General
=======

.. _cite:

How should I cite ASE?
----------------------

If you find ASE useful in your research please cite:

   | Ask Hjorth Larsen, Jens Jørgen Mortensen, Jakob Blomqvist,
   | Ivano E. Castelli, Rune Christensen, Marcin Dułak, Jesper Friis,
   | Michael N. Groves, Bjørk Hammer, Cory Hargus, Eric D. Hermes,
   | Paul C. Jennings, Peter Bjerre Jensen, James Kermode, John R. Kitchin,
   | Esben Leonhard Kolsbjerg, Joseph Kubal, Kristen Kaasbjerg,
   | Steen Lysgaard, Jón Bergmann Maronsson, Tristan Maxson, Thomas Olsen,
   | Lars Pastewka, Andrew Peterson, Carsten Rostgaard, Jakob Schiøtz,
   | Ole Schütt, Mikkel Strange, Kristian S. Thygesen, Tejs Vegge,
   | Lasse Vilhelmsen, Michael Walter, Zhenhua Zeng, Karsten Wedel Jacobsen
   | `The Atomic Simulation Environment—A Python library for working with atoms`__
   | J. Phys.: Condens. Matter Vol. **29** 273002, 2017

   __ https://doi.org/10.1088/1361-648X/aa680e

An older paper corresponding to an early version of ASE is:

   | S. R. Bahn and K. W. Jacobsen
   | `An object-oriented scripting interface to a legacy electronic structure code`__
   | Comput. Sci. Eng., Vol. **4**, 56-66, 2002

   __ http://dx.doi.org/10.1109/5992.998641

BibTex (:git:`doc/ASE.bib`):

.. literalinclude:: ASE.bib
