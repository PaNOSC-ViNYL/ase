.. module:: transport
   :synopsis: Electron transport

==================
Electron transport
==================

.. default-role:: math

The :mod:`transport` module of ASE assumes the generic setup of the system in
question sketched below:

. . . |setup| . . .

.. |setup| image:: transport_setup.png
   :align: middle

There is a central region (blue atoms plus the molecule) connected to
two semi-infinite leads constructed by infinitely repeated *principal
layers* (red atoms). The entire structure may be periodic in the
transverse direction, which can be effectively sampled using
**k**-points.

The system is described by a Hamiltonian matrix which must be
represented in terms of a localized basis set such that each element
of the Hamiltonian can be ascribed to either the left, central, or
right region, *or* the coupling between these.

The Hamiltonian can thus be decomposed as:

.. math::

    H = \begin{pmatrix}
      \ddots      & V_L         &             &             &     \\
      V_L^\dagger & H_L         & V_L         &             &     \\
                  & V_L^\dagger & H_C         & V_R         &     \\
                  &             & V_R^\dagger & H_R         & V_R \\
                  &             &             & V_R^\dagger & \ddots
    \end{pmatrix}

where `H_{L/R}` describes the left/right principal layer, and `H_C`
the central region. `V_{L/R}` is the coupling between principal
layers, *and* from the principal layers into the central region. The
central region must contain at least one principal layer on each side,
and more if the potential has not converged to its bulk value at this
size. The central region is assumed to be big enough that there is no
direct coupling between the two leads.

Having defined `H_{L/R}`, `V_{L/R}`, and `H_C`, the elastic
transmission function can be determined using the Nonequilibrium
Green Function (NEGF) method.  This is achieved by the class
:class:`~ase.transport.calculators.TransportCalculator` (in
ase.transport.calculators) which makes no requirement on the origin of
these five matrices.

For an example of how to use the :mod:`transport` module, see the GPAW
exercise on `electron transport`_

.. _electron transport: http://wiki.fysik.dtu.dk/gpaw/exercises/transport/transport.html

.. default-role::
