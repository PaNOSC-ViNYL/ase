.. _diffusion_tutorial:

=========================================
Diffusion of gold atom on Al(100) surface
=========================================

First, set up the initial and final states:

|initial| |final|

.. literalinclude:: diffusion1.py

.. note::  Notice how the tags are used to select the constrained atoms
   ...  XXX

Now, do the NEB calculation:

.. literalinclude:: diffusion2.py

|ts| |barrier|



.. |initial| image:: ../../_static/diffusion-I.png
.. |final| image:: ../../_static/diffusion-F.png
.. |ts| image:: ../../_static/diffusion-T.png
.. |barrier| image:: ../../_static/diffusion-barrier.png


Parallelizing over images
=========================

Instead of having one process do the calculations for all three
internal images in turn, it will be faster to have three processes do
one image each.  This can be done like this:

.. literalinclude:: diffusion3.py


