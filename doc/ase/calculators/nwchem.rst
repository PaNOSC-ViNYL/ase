.. module:: ase.calculators.nwchem

======
NWChem
======

`NWChem <http://www.nwchem-sw.org>`_ is a computational chemistry code
based on gaussian basis functions or plane-waves.


Setup
=====

.. highlight:: bash

You first need to install a working copy of NWChem for ASE to call;
follow the instructions on the `NWChem website <http://www.nwchem-sw.org>`_.

The default command that ASE wil use to start NWChem is
``nwchem PREFIX.nw > PREFIX.out``. You can change this command by setting the 
environment variable :envvar:`ASE_NWCHEM_COMMAND`. (For example, add a line
to your ``.bashrc`` with ``export ASE_NWCHEM_COMMAND="my new command"``.)

The default command will only allow you to run NWChem on a single core. To
run on multiple processors you will need to specify your MPI (or similar)
command, which typically requires telling the MPI command the number of tasks
to dedicate to the process. An example command to allow multiprocessing is
``mpirun -n $SLURM_NTASKS nwchem PREFIX.nw > PREFIX.out``, for the SLURM
queueing system. If you use a different queueing system replace ``$SLURM_NTASKS``
with the appropriate variable, such as ``$PBS_NP``.


Examples
========

Here is a command line example of how to calculate optimize the geometry of a
water molecule using PBE::

  $ asec H2O optimize -c nwchem -p xc=PBE 
  LBFGS:   0  16:17:29    -2064.914841       1.9673
  LBFGS:   1  16:17:31    -2064.963074       0.9482
  LBFGS:   2  16:17:32    -2064.976603       0.1425
  LBFGS:   3  16:17:33    -2064.977216       0.0823
  LBFGS:   4  16:17:35    -2064.977460       0.0010
  $ ase-gui H2O.traj@-1 -tg "a(1,0,2),d(0,1)"
  102.577881445 1.00806894632

.. highlight:: python

An example of creating an NWChem calculator in the python interface is::

  from ase.calculators.nwchem import NWChem
  
  calc = NWChem(label='calc/nwchem',
                maxiter=2000,
                xc='B3LYP',
                basis='6-31+G*')


Parameters
==========

The list of possible parameters and their defaults is shown below.

=================  =========  ========================= ============================
keyword            type       default value             description
=================  =========  ========================= ============================
``label``           ``str``    ``'nwchem'``             Label for saved files.
``xc``              ``str``    ``'LDA'``                Exchange-correlation functional.
``smearing``                   ``None``                 Smearing.
``charge``                     ``None``                 Charge
``task``            ``str``    ``'gradient'``           Task to perform. 'gradient' means force call.
``geometry``        ``str``    ``'nocenter noautosym'`` Geometry arguments. Note
                                                        NWChem centers the coordinates
                                                        by default.
``convergence``     ``dict``                            Convergence criteria.
``basis``           ``str``    ``'3-21G'``              Basic set.
``basispar``                   ``None``                 
``ecp``                        ``None``                 
``so``                         ``None``                 
``spinorbit``                  ``None``                 
``odft``                       ``None``                 
``raw``                        ``''``                   Raw text outside DFT block
                                                        control string.
=================  =========  ========================= ============================

See the source code link below for further details.

.. autoclass:: NWChem
