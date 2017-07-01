.. module:: ase.calculators.checkpoint

=============
Checkpointing
=============

Checkpointing adds restart and rollback capabilities to ASE scripts. It stores
the current state of the simulation (and its history) into an :mod:`ase.db`.
Something like what follows is found in many ASE scripts::

  if os.path.exists('atoms_after_relax.traj'):
      a = ase.io.read('atoms_after_relax.traj')
  else:
      ase.optimize.FIRE(a).run(fmax=0.01)
      ase.io.write('atoms_after_relax.traj')

The idea behind checkpointing is to replace this manual checkpointing
capability with a unified infrastructure.


Manual checkpointing
====================

The class :class:`Checkpoint` takes care of storing and retrieving
information from the database. This information *always* includes an
:class:`~ase.Atoms` object, and it can include attached information on
the internal state of the script.

.. autoclass:: ase.calculators.checkpoint.Checkpoint
    :members:
    :member-order: bysource

In order to use checkpointing, first create a Checkpoint object::

  from ase.calculators.checkpoint import Checkpoint
  CP = Checkpoint()

You can optionally choose a database filename. Default is ``checkpoints.db``.

Code blocks are wrapped into checkpointed regions::

  try:
      a = CP.load()
  except NoCheckpoint:
      ase.optimize.FIRE(a).run(fmax=0.01)
      CP.save(a)

The code block in the ``except`` statement is executed only if it has not yet
been executed in a previous run of the script. The :meth:`~Checkpoint.save`
statement stores all of its parameters to the database.

This is not yet much shorter than the above example. The checkpointing object
can, however, store arbitrary information along the :class:`~ase.Atoms`
object. Imagine we have computed elastic constants and don't want to recompute
them. We can then use::

  try:
      a, C = CP.load()
  except NoCheckpoint:
      C = fit_elastic_constants(a)
      CP.save(a, C)

Note that one parameter to :meth:`~Checkpoint.save` needs to be an
:class:`~ase.Atoms` object, the others can be arbitrary. The
:meth:`~Checkpoint.load` statement returns these parameters in the order they
were stored upon save. In the above example, the elastic constants are stored
attached to the atomic configuration. If the script is executed again after the
elastic constants have already been computed, it will skip that computation and
just use the stored value.

If the checkpointed region contains a single statement, such as the above,
there is a shorthand notation available::

  C = CP(fit_elastic_constants)(a)

Sometimes it is necessary to checkpoint an iterative loop. If the script
terminates within that loop, it is useful to resume calculation from the same
loop position::

  try:
      a, converged, tip_x, tip_y = CP.load()
  except NoCheckpoint:
      converged = False
      tip_x = tip_x0
      tip_y = tip_y0
  while not converged:
      ... do something to find better crack tip position ...
      converged = ...
      CP.flush(a, converged, tip_x, tip_y)

The above code block is an example of an iterative search for a crack tip
position. Note that the convergence criteria needs to be stored to the database
so the loop is not executed if convergence has been reached. The
:meth:`~Checkpoint.flush` statement overrides the last value stored to the
database.

As a rule :meth:`~Checkpoint.save` has to be used inside an
``except NoCheckpoint`` statement and :meth:`~Checkpoint.flush` outside.


Automatic checkpointing with the checkpoint calculator
======================================================

The :class:`CheckpointCalculator` is a shorthand for wrapping every single
energy/force evaluation in a checkpointed region. It wraps the actual
calculator.

.. autoclass:: ase.calculators.checkpoint.CheckpointCalculator
    :members:
    :member-order: bysource

Example usage::

  calc = ...
  cp_calc = CheckpointCalculator(calc)
  atoms.set_calculator(cp_calc)
  e = atoms.get_potential_energy()

The first call to :meth:`~ase.Atoms.get_potential_energy` does the actual
calculation, a rerun of the script will load energies and force from the
database. Note that this is useful for calculation where each energy evaluation
is slow (e.g. DFT), but not recommended for molecular dynamics with classical
potentials since every single time step will be dumped to the database. This
will generate huge files.
