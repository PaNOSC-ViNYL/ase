.. module:: optimize

======================
Structure optimization
======================

There are currently 4 different optimization algorithms available:
``QuasiNewton``, ``MDMin``, ``FIRE``, and ``GLBFGS``.

``MDMin`` and ``FIRE`` both use Newtonian dynamics with added
friction, to converge to an energy minimum, whereas ``QuasiNewton``
uses the forces of consecutive steps to dynamically update a Hessian
describing the curvature of the potential energy landscape. ``GLBFGS``
is an experimental optimizer designed for simultaneous update of the
images along a nudged elastic band trajectory.

All optimizer classes have the following structure::

  class Optimizer:
      def __init__(self, atoms, restart=None, logfile=None):
      def run(self, fmax=0.05, steps=100000000):


QuasiNewton
-----------

The ``QuasiNewton`` object is one of the minimizers in the ASE
package.  Let's try to use it to optimize the structure of a water
molecule.  We start with the experimental geometry::

  from ase import *
  d = 0.9575
  t = pi / 180 * 104.51
  water = Atoms('H2O',
                positions=[(d, 0, 0),
                           (d * cos(t), d * sin(t), 0),
                           (0, 0, 0)],
                calculator=EMT())
  dyn = QuasiNewton(water)
  dyn.run(fmax=0.05)
  QuasiNewton:   0        6.445801      51.6847
  QuasiNewton:   1        2.418583      27.2946
  QuasiNewton:   2        0.551767      12.1607
  QuasiNewton:   3       -0.039301       4.0520
  QuasiNewton:   4       -0.128045       0.8479
  QuasiNewton:   5       -0.132312       0.0397

When doing structure optimization, it is useful to write the
trajectory to a file, so that the progress of the optimization run can
be followed during or after the run::

  traj = PickleTrajectory('H2O.traj', 'w', water)
  dyn = QuasiNewton(water)
  dyn.attatch(traj.write)
  dyn.run(fmax=0.05)
  
Use the command ``ag H2O.traj`` to see what is going on (more here: ase.gui_).

The ``attach`` method takes an optional argument ``interval=n`` that can
be used to tell the structure optimizer object to write the
configuration to the trajectory file only every ``n`` steps.


Hessian ...

Restart ...


LBFGS
-----

...

FIRE
----

...

MDMin
-----

The MDmin algorithm is a modification of the usual velocity-Verlet
molecular dynamics algorithm.  Newtons second law is solved
numerically, but after each time step the dot product between the
forces and the momenta is checked.  If it is zero, the system has just
passed through a (local) minimum in the potential energy, the kinetic
energy is large and about to decrease again.  At this point, the
momentum is set to zero.  Unlike a "real" molecular dynamics, the
masses of the atoms are not used, instead all masses are set to one.

The MDmin algorithm exists in two flavors, one where each atom is
tested and stopped individually (QuickMinAtomByAtom in the old ASE),
and one where all coordinates are treated as one long vector, and all
momenta are set to zero if the dotproduct between the momentum vector
and force vector (both of length 3N) is zero (QuickMinAllCoordinates
in the old ASE).  This module implements the latter version.

Although the algorithm is primitive, it performs very well because it
takes advantage of the physics of the problem.  Once the system is so
near the minimum that the potential energy surface is approximately
quadratic it becomes advantageous to switch to a minimization method
with quadratic convergence, such as `Conjugate Gradient` or `Quasi
Newton`.
