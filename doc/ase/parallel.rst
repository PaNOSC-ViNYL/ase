.. module:: ase.parallel

=====================
Parallel calculations
=====================

ASE will automatically run in parallel, if it can import an MPI communicator from any of the supported libraries. ASE will attepmt to import communicators from external libraries in following order: GPAW, Asap, Scientific MPI and MPI4PY.

If a parallel library is found, the ase.io.read function will always read only on master (of the MPI world object) and broadcast the atoms to all other cores. Therefore, always when using ase.io.read, all cores must read the same atoms in same order, for example in the case of NEB calculation.

If one requires an individual core/cores to read a particular file, please use Trajectory:

>>> from ase.io import Trajectory
>>> from gpaw.mpi import world
>>> atoms = Trajectory('myfile_%d.traj' % world.rank)[-1]

.. autofunction:: paropen
.. autofunction:: parprint
.. autofunction:: broadcast
.. autofunction:: parallel_function
.. autofunction:: parallel_generator
