# VelocityDistributions.py -- set up a velocity distribution

"""Module for setting up e.g. Maxwell-Boltzmann velocity distributions.

Currently, only one function is defined, MaxwellBoltzmannDistribution,
which sets the momenta of a list of atoms according to a
Maxwell-Boltzmann distribution at a given temperature.
"""

import numpy as np

def _maxwellboltzmanndistribution(masses, temp):
    xi = np.random.standard_normal((len(masses),3))
    momenta = xi * np.sqrt(masses * temp)[:,np.newaxis]
    return momenta

def MaxwellBoltzmannDistribution(atoms, temp):
    """Sets the momenta to a Maxwell-Boltzmann distribution."""
    momenta = _maxwellboltzmanndistribution(atoms.get_masses(), temp)
    atoms.set_momenta(momenta)


