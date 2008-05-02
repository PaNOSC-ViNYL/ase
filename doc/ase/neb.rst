===================
Nudged elastic band
===================

.. module:: neb
   :synopsis: Nudged Elastic Band method.

The Nudged Elastic Band method is a technique for finding transition paths
(and corresponding energy barriers) between given initial and final states.
The method involves constructing a "chain" of "replicas" or "images" of the
system and relaxing them in a certain way.

Relevant literature References:


1. H. Jonsson, G. Mills, and K. W. Jacobsen, in 'Classical and Quantum
   Dynamics in Condensed Phase Systems', edited by B. J. Berne,
   G. Cicotti, and D. F. Coker, World Scientific, 1998 [standard
   formulation]

2. 'Improved Tangent Estimate in the NEB method for Finding Minimum
   Energy Paths and Saddle Points', Graeme Henkelman and Hannes
   Jonsson, J. Chem. Phys. 113, 9978 (2000) [improved tangent
   estimates]

3. 'A Climbing-Image NEB Method for Finding Saddle Points and Minimum
   Energy Paths', Graeme Henkelman, Blas P. Uberuaga and Hannes
   Jonsson, J. Chem. Phys. 113, 9901 (2000)


The NEB class
=============

This module defines one class:

.. class:: NEB(images, k=0.1, climb=False)

Example of use::

  # Read initial and final states:
  initial = read('A.xyz')
  final = read('B.xyz')
  # Make a band consisting of 5 images:
  images = [initial]
  images += [initial.copy() for i in range(3)]
  images += [final]
  neb = NEB(images)
  # Interpolate linearly the potisions of the three middle images:
  neb.interpolate()
  # Set calculators:
  for image in images:
      image.set_calculator(MyCalculator(...))
  # Optimize:
  optimizer = QuasiNewton(neb)
  optimizer.run(fmax=0.04)

Notice the use of the :meth:`~NEB.interpolate` method to get a good
initial guess for the path from A to B.

.. method:: NEB.interpolate()

   Interpolate path linearly from initial to final state.

.. seealso::

   :mod:`optimize`:
        Information about energy minimization (optimization).

   :mod:`calculators`:
        How to use calculators.

   :ref:`tutorials`:

        * :ref:`neb1`
        * :ref:`neb2`



Trajectories
============

XXX


Climbing image
==============

The "climbing image" variation involves designating a specific image to behave
differently to the rest of the chain: it feels no spring forces, and the
component of the potential force parallel to the chain is reversed, such that
it moves towards the saddle point. This depends on the adjacent images
providing a reasonably good approximation of the correct tangent at the
location of the climbing image; thus in general the climbing image is not
turned on until some iterations have been run without it (generally 20% to 50%
of the total number of iterations).

To use the climbing image NEB method, instantiate the NEB object like this::

  neb = NEB(images, climb=True)
