.. module:: infrared


====================
Infrared intensities
====================

:class:`~ase.infrared.InfraRed` is an extension of
:class:`~ase.vibrations.Vibrations`, in addition to the
vibrational modes, also the infrared intensities of the modes
are calculated for an :class:`~ase.atoms.Atoms` object.

.. autoclass:: ase.infrared.InfraRed
   :members:

Efficient usage with VASP
=========================

Using VASP as calculator, to optimize the efficiency of the IR object the user 
should generate highly converged :file:`WAVECAR` and :file:`CHGCAR` files
for the optimized ionic structure of the system of interest. VASP should 
be set up to read these files at each displacement, but not re-write 
them at the end of a displacement, i.e. the :class:`~ase.calculators.vasp.Vasp`
calculator object should have the additional settings::

    >>> calc.set(istart = 1,
             	 icharg = 1,
             	 lwave = False,
             	 lcharg = False)

This way, each displacement will start with the wavefunctions from
the optimized structure, and the charge density extrapolated from
the optimized positions to the new positions of the ions, including
the displacement.