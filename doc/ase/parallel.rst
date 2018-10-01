.. module:: ase.parallel

=====================
Parallel calculations
=====================

ASE will automatically run in parallel, if it can import an MPI communicator
from any of the supported libraries. ASE will attempt to import communicators
from these external libraries: GPAW, Asap, Scientific MPI and
MPI4PY.

If a parallel library is found, the :func:`ase.io.read` function will always
read only on master (of the MPI world object) and broadcast the atoms to all
other cores. Therefore, always when using :func:`ase.io.read`, all cores must
read the same atoms in same order, for example in the case of a NEB
calculation.

If one requires an individual core/cores to read a particular file, please
use :func:`~ase.io.Trajectory`:

>>> from ase.io import Trajectory
>>> from ase.parallel import world
>>> atoms = Trajectory('myfile_{}.traj'.format(world.rank))[-1]

.. autofunction:: paropen
.. autofunction:: parprint
.. autofunction:: broadcast
.. autofunction:: parallel_function
.. autofunction:: parallel_generator
