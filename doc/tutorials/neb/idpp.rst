.. _idpp_tutorial:

==============================================================================
Image Dependent Pair Potential for improved interpolation of NEB initial guess
==============================================================================

Reference: S. Smidstrup, A. Pedersen, K. Stokbro and H. Jonsson,
Improved initial guess for minimum energy path calculations,
J. Chem. Phys. 140, 214106 (2014).

Use of the NEB method is dependent upon generating an initial guess 
for the images lying between the initial and final states. The most 
simple approach is to use linear interpolation of the 
atomic coordinates. However, this can be problematic as the quality 
of the interpolated path can ofter be far from the real one. The implication being
that a lot of time is spent in the NEB routine optimising the shape of
the path, before the transition state is homed-in upon. 

The image dependent pair potential is a method that has been 
developed to provide an improvement to the initial guess for the NEB path. 
The IDPP method uses the bond distance between the atoms involved in 
the transition state to create target structures for the images, rather
than interpolating the atomic positions. By defining an objective function in terms
of the distances between atoms, the NEB algorithm is used with this 
image dependent pair potential (IDPP) to create the initial guess for the
full NEB calculation.

Note: The examples below utilise the EMT calculator for illustrative purposes, the
results should not be over interpreted. 

Example 1: Ethane 
=================

This example illustrates the use of the IDPP interpolation scheme to 
generate an initial guess for rotation of a methyl group around the CC bond. 

Using the standard linear interpolation approach, as in the following example, we can see
that 46 iterations are required to find the transition state.

.. literalinclude:: idpp1.py

However if we modify our script slightly and use the IDPP method to find the initial
guess, we can see that the number of iterations required to find the transition
state is reduced to 14. 

.. literalinclude:: idpp2.py

Clearly, if one was using a full DFT calculator one can
potentially gain a significant time improvement.

Example 2: N diffusion over a step edge
=======================================

Often we are interested in generating an initial guess for a surface reaction. 
This example illustrates how we can optimise our initial and final state structures
before using the IDPP interpolation to generate our initial guess 
for the NEB calculation:

.. literalinclude:: idpp3.py

To again illustrate the potential speedup, the following script which
uses the linear interpolation takes 23 iterations to find a transition
state, compared to 12 using the IDPP interpolation. 

.. literalinclude:: idpp4.py


