# VelocityDistributions.py -- set up a velocity distribution

"""Module for setting up e.g. Maxwell-Boltzmann velocity distributions.

Currently, only one function is defined, MaxwellBoltzmannDistribution,
which sets the momenta of a list of atoms according to a
Maxwell-Boltzmann distribution at a given temperature.
"""

import numpy as np

# For parallel GPAW simulations, the random velocities should be distributed.
if '_gpaw' in sys.modules:
    # http://wiki.fysik.dtu.dk/gpaw
    from gpaw.mpi import world as gpaw_world
else:
    gpaw_world = None

def _maxwellboltzmanndistribution(masses, temp):
    xi = np.random.standard_normal((len(masses),3))
    if gpaw_world is not None:
        gpaw_world.broadcast(xi, 0)
    momenta = xi * np.sqrt(masses * temp)[:,np.newaxis]
    return momenta

def MaxwellBoltzmannDistribution(atoms, temp):
    """Sets the momenta to a Maxwell-Boltzmann distribution."""
    momenta = _maxwellboltzmanndistribution(atoms.get_masses(), temp)
    atoms.set_momenta(momenta)

def Stationary(atoms):
    "Sets the center-of-mass momentum to zero."
    p = atoms.get_momenta()
    p0 = np.sum(p, 0)
    # We should add a constant velocity, not momentum, to the atoms
    m = atoms.get_masses()
    mtot = np.sum(m)
    v0 = p0/mtot
    p -= v0*m[:,np.newaxis]
    atoms.set_momenta(p)
    
