.. module:: ase.calculators.vasp

====
VASP
====

Introduction
============

VASP_ is a density-functional theory code using pseudopotentials or 
the projector-augmented wave method and a plane wave basis set. This 
interface makes it possible to use VASP_ as a calculator in ASE, and 
also to use ASE as a post-processor for an already performed VASP_
calculation.


.. _VASP: http://cms.mpi.univie.ac.at/vasp/



Environment variables
=====================

You need to write a script called :file:`run_vasp.py` containing
something like this::

  import os
  exitcode = os.system('vasp')

The environment variable :envvar:`VASP_SCRIPT` must point to that file.

A directory containing the pseudopotential directories :file:`potpaw` 
(LDA XC) :file:`potpaw_GGA` (PW91 XC) and :file:`potpaw_PBE` (PBE XC)
is also needed, and it is to be put in the environment variable
:envvar:`VASP_PP_PATH`.

Set both environment variables in your shell configuration file:

.. highlight:: bash
 
::

  $ export VASP_SCRIPT=$HOME/vasp/run_vasp.py
  $ export VASP_PP_PATH=$HOME/vasp/mypps

.. highlight:: python



VASP Calculator
=============== 

The default setting used by the VASP interface is

.. class:: Vasp(restart=None, xc='PW91', setups=None, kpts=(1,1,1), gamma=None)

Below follows a list with a selection of parameters

==============  =========  ==============  ============================
keyword         type       default value   description
==============  =========  ==============  ============================
``restart``     ``bool``   None            Restart old calculation or
                                           use ASE for post-processing
``xc``          ``str``    'PW91'          XC-functional. Defaults to
                                           None if ``gga`` set explicitly.
``setups``      ``str``    None            Additional setup option
``pp``          ``str``    Set by ``xc``   Pseudopotential (POTCAR) set
                           or ``gga``      used (LDA, PW91 or PBE).
``kpts``        *seq*      `\Gamma`-point  **k**-point sampling
``gamma``       ``bool``   None            `\Gamma`-point centered 
                                           **k**-point sampling
``reciprocal``  ``bool``   None            Use reciprocal units if
                                           **k**-points are specified 
                                           explicitly
``prec``        ``str``                    Accuracy of calculation
``encut``       ``float``                  Kinetic energy cutoff
``ediff``       ``float``                  Convergence break condition
                                           for SC-loop.
``nbands``      ``int``                    Number of bands
``algo``        ``str``                    Electronic minimization 
                                           algorithm
``ismear``      ``int``                    Type of smearing
``sigma``       ``float``                  Width of smearing
``nelm``        ``int``                    Maximum number of
                                           SC-iterations
==============  =========  ==============  ============================

*seq*: A sequence of three ``int``'s.

For parameters in the list without default value given, VASP will set 
the default value. Most of the parameters used in the VASP :file:`INCAR` file 
are allowed keywords. See the official `VASP manual`_ for more details.

.. _VASP manual: http://cms.mpi.univie.ac.at/vasp/vasp/vasp.html


.. note::

   Parameters can be changed after the calculator has been constructed
   by using the :meth:`~ase.calculators.vasp.Vasp.set` method:

   >>> calc.set(prec='Accurate', ediff=1E-5)

   This would set the precision to Accurate and the break condition
   for the electronic SC-loop to ``1E-5`` eV.

Exchange-correlation functionals
================================

The ``xc`` parameter is used to define a "recipe" of other parameters
including the pseudopotential set ``pp``.  It is possible to override
any parameters set with ``xc`` by setting them explicitly. For
example, the screening parameter of a HSE calculation might be
modified with

   >>> calc = ase.calculators.vasp.Vasp(xc='hse06', hfscreen=0.4)

The default pseudopotential set is potpaw_PBE unless ``xc`` or ``pp``
is set to ``pw91`` or ``lda``.

==========================  =====================================
``xc`` value                Parameters set
==========================  =====================================
lda, pbe, pw91              ``pp`` (``gga`` set implicity in POTCAR)
pbesol, revpbe, rpbe, am05  ``gga``
tpss, revtpss, m06l         ``metagga``
vdw-df, optpbe-vdw          ``gga``, ``luse_vdw``, ``aggac``
optb88-vdw, obptb86b-vdw    ``gga``, ``luse_vdw``, ``aggac``,
                            ``param1``, ``param2``
beef-vdw                    ``gga``, ``luse_vdw``, ``zab_vdw``
vdw-df2                     ``gga``, ``luse_vdw``, ``aggac``,
                            ``zab_vdw``
hf                          ``lhfcalc``, ``aexx``, ``aldac``,
                            ``aggac``
pbe0                        ``gga``, ``lhfcalc``
b3lyp                       ``gga``, ``lhfcalc``, ``aexx``, ``aggax``,
                            ``aggac``, ``aldac``
hse03, hse06, hsesol        ``gga``, ``lhfcalc``, ``hfscreen``
==========================  =====================================

It is possible for the user to temporarily add their own ``xc``
recipes without modifying ASE, by updating a dictionary. For example,
to implement a hybrid PW91 calculation:

.. code-block:: python

   from ase.calculators.vasp import Vasp
   Vasp.xc_defaults['pw91_0'] = {'gga': '91', 'lhfcalc': True}

   calc = Vasp(xc='PW91_0')

Note that the dictionary keys must be *lower case*, while the ``xc``
parameter is case-insensitive when used.


Setups
======

For many elements, VASP is distributed with a choice of
pseudopotential setups. These may be hard/soft variants of the
pseudopotential or include additional valence electrons. While the
Vasp calculator will default to the pseudopotential folders with the
same name as the element, alternative setups may be selected
with the `setups` dictionary.

To use an alternative setup for all instances of an element, simply
provide the characters which need to be added, e.g.

.. code-block:: python
                
   calc = Vasp(xc='PBE', setups={'Li': '_sv'})

will use the ``Li_sv`` all-electron pseudopotential for all Li atoms.
To apply special setups to individual atoms, identify them by their
zero-indexed number in the atom list and use the full setup name. For
example,

.. code-block:: python
                
   calc= Vasp(xc='PBE', setups={3: 'Ga_d'})

will treat the Ga atom in position 3 (i.e. the fourth atom) of the
atoms object as special, with an additional 10 d-block valence
electrons, while other Ga atoms use the default 3-electron setup and
other elements use their own default setups. The positional index may
be quoted as a string (e.g. ``{'3': 'Ga_d'}``).

Spin-polarized calculation
==========================

If the atoms object has non-zero magnetic moments, a spin-polarized
calculation will be performed by default.

Here follows an example how to calculate the total magnetic moment of
a sodium chloride molecule.
  
.. literalinclude:: NaCl.py

In this example the initial magnetic moments are assigned to the atoms
when defining the Atoms object. The calculator will detect that at least
one of the atoms has a non-zero magnetic moment and a spin-polarized
calculation will automatically be performed. The ASE generated :file:`INCAR`
file will look like:

.. literalinclude:: INCAR_NaCl


.. note:: 
   
   It is also possible to manually tell the calculator to perform a 
   spin-polarized calculation:

   >>> calc.set(ispin=2)

   This can be useful for continuation jobs, where the initial magnetic 
   moment is read from the WAVECAR file.


Restart old calculation
=======================

To continue an old calculation which has been performed without the interface
use the ``restart`` parameter when constructing the calculator

>>> calc = Vasp(restart=True)

Then the calculator will read atomic positions from the :file:`CONTCAR` file,
physical quantities from the :file:`OUTCAR` file, **k**-points from the
:file:`KPOINTS` file and parameters from the :file:`INCAR` file. 

.. note::

   Only Monkhorst-Pack and \Gamma-centered **k**-point sampling are supported 
   for restart at the moment. Some :file:`INCAR` parameters may not be
   implemented for restart yet. Please report any problems to the ASE mailing
   list.

The ``restart`` parameter can be used , as the name suggest to continue a job from where a
previous calculation finished. Furthermore, it can be used to extract data from
an already performed calculation. For example, to get the total potential energy
of the sodium chloride molecule in the previous section, without performing any additional
calculations, in the directory of the previous calculation do:

>>> calc = Vasp(restart=True)
>>> atoms = calc.get_atoms()
>>> atoms.get_potential_energy()
-4.7386889999999999
