# VelocityDistributions.py -- set up a velocity distribution

"""Module for setting up e.g. Maxwell-Boltzmann velocity distributions.

Currently, only one function is defined, MaxwellBoltzmannDistribution,
which sets the momenta of a list of atoms according to a
Maxwell-Boltzmann distribution at a given temperature.
"""

import sys
import numpy as np
from ase.parallel import world


def _maxwellboltzmanndistribution(masses, temp, communicator=world):
    # For parallel GPAW simulations, the random velocities should be
    # distributed.  Uses gpaw world communicator as default, but allow
    # option of specifying other communicator (for ensemble runs)
    xi = np.random.standard_normal((len(masses), 3))
    communicator.broadcast(xi, 0)
    momenta = xi * np.sqrt(masses * temp)[:, np.newaxis]
    return momenta


def MaxwellBoltzmannDistribution(atoms, temp, communicator=world):
    """Sets the momenta to a Maxwell-Boltzmann distribution."""
    momenta = _maxwellboltzmanndistribution(atoms.get_masses(), temp,
                                            communicator)
    atoms.set_momenta(momenta)


def Stationary(atoms):
    "Sets the center-of-mass momentum to zero."
    p = atoms.get_momenta()
    p0 = np.sum(p, 0)
    # We should add a constant velocity, not momentum, to the atoms
    m = atoms.get_masses()
    mtot = np.sum(m)
    v0 = p0 / mtot
    p -= v0 * m[:, np.newaxis]
    atoms.set_momenta(p)
