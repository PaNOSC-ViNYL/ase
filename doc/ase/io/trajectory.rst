.. module:: ase.io.trajectory
   :synopsis: Trajectory input-output module

================
Trajectory files
================

.. contents::
   
The :mod:`ase.io.trajectory` module defines Trajectory objects, that is
objects storing the temporal evolution of a simulation or the path
taken during an optimization.  A Trajectory file
contains one or more :class:`~ase.Atoms` objects, usually to be
interpreted as a time series, although that is not a requirement.

The main Trajectory object writes in a file format, which is compatible
across Python version.

The :mod:`ase.io.trajectory` additionally defines two specialized kinds of
Trajectory files, the PickleTrajectory and the BundleTrajectory.

PickleTrajectory is the old (pre 2015) Trajectory format, its use is no
longer recommended as compatibility between Python versions (and to a
lesser degree between ASE vesions) cannot be guaranteed.  You *must*
:ref:`convert your old PickleTrajectory files <convert>` as soon
as possible.

BundleTrajectory is only intended for large molecular dynamics
simulations (large meaning millions of atoms).

Typically, trajectories are used to store different configurations of
the same system (i.e. the same atoms).  If you need to store
configurations of different systems, the :mod:`ASE Database module
<ase.db>` may be more appropriate.


Trajectory
==========

The Trajectory function returns a Trajectory reading or writing
object, depending on the mode.

.. autofunction:: ase.io.Trajectory

The function returns a TrajectoryReader or a TrajectoryWriter object.

Reading a trajectory file is done by indexing the TrajectoryReader
object, i.e. traj[0] reads the first configuration, traj[-1] reads the
last, etc.

Writing a trajectory file is done by calling the ``write`` method.  If no
atoms object was given when creating the object, it must be given as
an argument to the ``write`` method.


Examples
--------

Reading a configuration::

    from ase.io.trajectory import Trajectory
    traj = Trajectory('example.traj')
    atoms = traj[-1]

Reading all configurations::

    traj = Trajectory('example.traj')
    for atoms in traj:
        # Analyze atoms

Writing every 100th time step in a molecular dynamics simulation::

    # dyn is the dynamics (e.g. VelocityVerlet, Langevin or similar)
    traj = Trajectory('example.traj', 'w', atoms)
    dyn.attach(traj.write, interval=100)
    dyn.run(10000)
    traj.close()

    
.. _new trajectory:
    
The TrajectoryReader and TrajectoryWriter objects
-------------------------------------------------

Usually, you only need the interface given above, but the reader and
writer have a few additional methods, that can be useful.

.. autoclass:: ase.io.trajectory.TrajectoryReader
   :members:

Note that there is apparently no methods for reading the trajectory.
Reading is instead done by indexing the trajectory, or by iterating
over the trajectory: ``traj[0]`` and ``traj[-1]`` return the first and
last :class:`~ase.Atoms` object in the trajectory.

.. autoclass:: ase.io.trajectory.TrajectoryWriter
   :members:


.. _old trajectory:
      
PickleTrajectory
================

The *obsolete* PickleTrajectory uses the same object for reading and writing.

**WARNING 1:** If your Atoms objects contains constraints, the
constraint object is pickled and stored in the file.  Unfortunately,
this means that if the object definition in ASE changes, you cannot
read the trajectory file.  In the new
Trajectory format the contraint is stored in an
implementation-independent format.

**WARNING 2:** It is possible to write a malicious pickle file (and
thus a malicious PickleTrajectory) that executes arbitrary code when
reading the file.  The new Trajectory format cannot contain code.

For the reasons above, version 3.10 of ASE will not be able to read and write
PickleTrajectory files, and you need to :ref:`convert existing files <convert>`
to the new format.

      
.. _convert:

Converting old PickleTrajectory files to new Trajectory files
-------------------------------------------------------------

Please convert your old PickleTrajectory files before it is too late::

    $ python -m ase.io.trajectory file1.traj [file2.traj ...]

this will convert one or more files.  The original files are kept with
extension ``.traj.old``

You can identify old trajectory files like this::

    $ python -m ase.io.formats hmmm.traj
    hmmm.traj: Old ASE pickle trajectory (trj+)
    $ python -m ase.io.trajectory hmmm.traj  # convert
    $ python -m ase.io.formats hmmm.traj hmmm.traj.old
    hmmm.traj:     ASE trajectory (traj+)
    hmmm.traj.old: Old ASE pickle trajectory (trj+)


BundleTrajectory
================

The BundleTrajectory has the interface

.. autoclass:: ase.io.bundletrajectory.BundleTrajectory
   :members:

      
See also
========

* The function :func:`ase.io.write` can write a single
  :class:`~ase.Atoms` object to a Trajectory file.

* The function :func:`ase.io.read` can read an :class:`~ase.Atoms`
  object from a Trajectory file, per default it reads the last one.

* The database modue :mod:`ase.db`.

