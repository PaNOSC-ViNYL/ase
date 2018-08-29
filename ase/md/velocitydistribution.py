# encoding: utf-8
# VelocityDistributions.py -- set up a velocity distribution

"""Module for setting up velocity distributions such as Maxwell–Boltzmann.

Currently, only a few functions are defined, such as
MaxwellBoltzmannDistribution, which sets the momenta of a list of
atoms according to a Maxwell-Boltzmann distribution at a given
temperature.

"""

import numpy as np
from ase.parallel import world


def _maxwellboltzmanndistribution(masses, temp, communicator=world,
                                  rng=np.random):
    # For parallel GPAW simulations, the random velocities should be
    # distributed.  Uses gpaw world communicator as default, but allow
    # option of specifying other communicator (for ensemble runs)
    xi = rng.standard_normal((len(masses), 3))
    communicator.broadcast(xi, 0)
    momenta = xi * np.sqrt(masses * temp)[:, np.newaxis]
    return momenta


def MaxwellBoltzmannDistribution(atoms, temp, communicator=world,
                                 force_temp=False, rng=np.random):
    """Sets the momenta to a Maxwell-Boltzmann distribution. temp should be
    fed in energy units; i.e., for 300 K use temp=300.*units.kB. If
    force_temp is set to True, it scales the random momenta such that the
    temperature request is precise.
    """
    momenta = _maxwellboltzmanndistribution(atoms.get_masses(), temp,
                                            communicator, rng)
    atoms.set_momenta(momenta)
    if force_temp:
        temp0 = atoms.get_kinetic_energy() / len(atoms) / 1.5
        gamma = temp / temp0
        atoms.set_momenta(atoms.get_momenta() * np.sqrt(gamma))


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


def ZeroRotation(atoms):
    "Sets the total angular momentum to zero by counteracting rigid rotations."
    # Find the principal moments of inertia and principal axes basis vectors
    Ip, basis = atoms.get_moments_of_inertia(vectors=True)
    # Calculate the total angular momentum and transform to principal basis
    Lp = np.dot(basis, atoms.get_angular_momentum())
    # Calculate the rotation velocity vector in the principal basis, avoiding
    # zero division, and transform it back to the cartesian coordinate system
    omega = np.dot(np.linalg.inv(basis), np.select([Ip > 0], [Lp / Ip]))
    # We subtract a rigid rotation corresponding to this rotation vector
    com = atoms.get_center_of_mass()
    positions = atoms.get_positions()
    positions -= com  # translate center of mass to origin
    velocities = atoms.get_velocities()
    atoms.set_velocities(velocities - np.cross(omega, positions))


def PhononHarmonics(atoms, temp, force_constants, rng=np.random,
                    failfast=True):
    r"""Excite phonon modes to specified temperature.

    This excites all phonon modes randomly so that each contributes,
    on average, equally to the given temperature.  Both potential
    energy and kinetic energy will be consistent with the phononic
    vibrations characteristic of the specified temperature.

    In other words the system will be equilibrated for an MD run at
    that temperature.

    force_constants should be the matrix as force constants, e.g.,
    as computed by the ase.phonons module.

    Let X_ai be the phonon modes indexed by atom and mode, w_i the
    phonon frequencies, and let 0 < Q_i <= 1 and 0 <= R_i < 1 be
    uniformly random numbers.  Then

                   1/2
      _     / k T \     ---  1  _             1/2
      R  += | --- |      >  --- X   (-2 ln Q )    cos (2 pi R )
       a    \  m  /     ---  w   ai         i                i
                a        i    i


                   1/2
      _     / k T \     --- _            1/2
      v   = | --- |      >  X  (-2 ln Q )    sin (2 pi R )
       a    \  m  /     ---  ai        i                i
                a        i

    See e.g.:

      http://ollehellman.github.io/program/canonical_configuration.html

    for a derivation.

    If failfast is true, a ValueError will be raised if the
    eigenvalues of the dynamical matrix look suspicious (they should
    be three zeros for translational modes followed by any number of
    strictly positive eigenvalues)."""
    mass_a = atoms.get_masses()

    rminv = (mass_a**-0.5).repeat(3)
    dynamical_matrix = force_constants.copy()
    dynamical_matrix *= rminv[:, None]
    dynamical_matrix *= rminv[None, :]

    w2_s, X_is = np.linalg.eigh(dynamical_matrix)

    if failfast:
        zeros = w2_s[:3]
        worst_zero = np.abs(zeros).max()
        if worst_zero > 1e-3:
            raise ValueError('Translational modes have suspiciously large '
                             'energies; should be close to zero: {}'
                             .format(w2_s[:3]))

        w2min = w2_s[3:].min()
        if w2min < 0:
            raise ValueError('Dynamical matrix has negative eigenvalues '
                             'such as {}'.format(w2min))

    # First three modes are translational so ignore:
    nw = len(w2_s) - 3
    w_s = np.sqrt(w2_s[3:])
    X_acs = X_is[:, 3:].reshape(len(atoms), 3, nw)

    # We need 0 < P <= 1.0 and not 0 0 <= P < 1.0 for the logarithm
    # to avoid (highly improbable) NaN:
    A_s = np.sqrt(-2.0 * np.log(1.0 - rng.rand(nw)))
    phi_s = 2.0 * np.pi * rng.rand(nw)

    d_ac = (A_s * np.sin(phi_s) / w_s * X_acs).sum(axis=2)
    v_ac = (A_s * np.cos(phi_s) * X_acs).sum(axis=2)

    sqrtkTm_a = np.sqrt(temp / mass_a)
    for quantity in [d_ac, v_ac]:
        quantity *= sqrtkTm_a[:, None]

    atoms.positions += d_ac
    atoms.set_velocities(v_ac)
