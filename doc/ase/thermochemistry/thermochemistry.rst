.. module:: ase.thermochemistry
   :synopsis: Thermochemistry module

===============
Thermochemistry
===============

ASE contains a :mod:`ase.thermochemistry` module that lets the user derive
commonly desired thermodynamic quantities of molecules and crystalline solids
from ASE output and some user-specified parameters. Four cases are currently
handled by this module: the ideal-gas limit (in which translational and
rotational degrees of freedom are taken into account), the hindered
translator / hindered rotor model (used for adsorbates, in which two degrees
of freedom are translational, one is rotational, and the remaining 3N-3 are
vibrational), the harmonic limit (generally used for adsorbates, in which all
degrees of freedom are treated harmonically), and a crystalline solid model
(in which a lattice of N atoms is treated as a system of 3N independent
harmonic oscillators). The first three cases rely on good vibrational energies
being fed to the calculators, which can be calculated with the
:mod:`ase.vibrations` module. Likewise, the crystalline solid model depends on
an accurate phonon density of states; this is readily calculated using the
:mod:`ase.phonons` module.


Ideal-gas limit
===============

The thermodynamic quantities of ideal gases are calculated by assuming that
all spatial degrees of freedom are independent and separable into
translational, rotational, and vibrational degrees of freedom. The
:class:`~ase.thermochemistry.IdealGasThermo` class supports calculation of
enthalpy (:math:`H`), entropy (:math:`S`), and Gibbs free energy (:math:`G`),
and has the interface listed below.

.. autoclass:: IdealGasThermo
   :members:


Example
-------

The :class:`IdealGasThermo` class would generally be called after an energy
optimization and a vibrational analysis. The user needs to supply certain
parameters if the entropy or free energy are desired, such as the geometry
and symmetry number. An example on the nitrogen molecule is:

.. literalinclude:: nitrogen.py


This will give the thermodynamic summary output:

.. literalinclude:: nitrogen.txt


Hindered translator / hindered rotor model
==========================================

The hindered translator / hindered rotor model bridges the gap between the
2D gas (i.e. free translator / free rotor) and the 2D lattice gas (i.e.
harmonic oscillator). For an adsorbate containing N atoms, two degrees of
freedom are treated as hindered translations in the two directions parallel to
the surface, one degree of freedom is treated as a hindered rotation about the
axis perpendicular to the surface, and the remaining 3N-3 degrees of freedom
are treated as vibrations. The :class:`HinderedThermo` class supports the
calculation of internal energy, entropy, free energy, and zero point energy
(included in the internal energy). All of the thermodynamic properties
calculated here are at the standard state surface concentration (defined here
such that a 2D ideal gas at that concentration has 2/3 the translational
entropy of a 3D ideal gas at 1 bar pressure, so that :math:`\theta^0` = 0.012
at 298 K for a surface with `10^{15}` sites/cm\ :sup:`2`). This class
returns the Helmholtz free energy; if the user assumes that the pV term (in G =
U + pV - TS) is zero then this free energy can also be interpreted as the Gibbs
free energy. This class depends on the user defined translation barrier
(trans_barrier_energy) and rotational barrier (rot_barrier_energy) for the
adsorbate to move on the surface in order to calculate the translational and
rotational degrees of freedom. To calculate the vibrational degrees of freedom,
all 3N vibrational energies must be supplied in the vib_energies list and the
3N-3 largest vibrational energies are used to calculate the vibrational
contribution; this is a list as can be generated with the .get_energies()
method of :class:`ase.vibrations.Vibrations`. The class :class:`HinderedThermo`
has the interface described below.

.. autoclass:: HinderedThermo
   :members:


Example
-------

The :class:`HinderedThermo` class would generally be called after an energy
optimization and a vibrational analysis. The user needs to supply certain
parameters, such as the vibrational energies, translational energy barrier,
rotational energy barrier, surface site density, number of equivalent minima
in a full rotation, and the number of symmetric arms of the adsorbate as it
rotates on the surface. The user also needs to supply either the mass of the
adsorbate and the reduced moment of inertia of the adsorbate as it rotates on
the surface or the user can supply the atoms object from which the mass and an
approximate reduced moment of inertia may be determined. An example for ethane
on a platinum (111) surface is:

.. literalinclude:: ethane.py


This will give the thermodynamic summary output:

.. literalinclude:: ethane.txt


Harmonic limit
==============

In the harmonic limit, all degrees of freedom are treated harmonically. The
:class:`HarmonicThermo` class supports the calculation of internal energy,
entropy, and free energy. This class returns the Helmholtz free energy; if
the user assumes the pV term (in H = U + pV) is zero this can also be
interpreted as the Gibbs free energy. This class uses all of the energies
given to it in the vib_energies list; this is a list as can be generated
with the .get_energies() method of :class:`ase.vibrations.Vibrations`, but
the user should take care that all of these energies are real
(non-imaginary). The class :class:`HarmonicThermo` has the interface
described below.

.. autoclass:: HarmonicThermo
   :members:


Crystals
========

In this model a crystalline solid is treated as a periodic system of
independent harmonic oscillators. The :class:`CrystalThermo` class supports
the calculation of internal energy (:math:`U`), entropy (:math:`S`) and
Helmholtz free energy (:math:`F`), and has the interface listed below.

.. autoclass:: CrystalThermo
   :members:


Example
-------

The :class:`CrystalThermo` class will generally be called after an energy
optimization and a phonon vibrational analysis of the crystal. An example for
bulk gold is:

.. literalinclude:: gold.py

This will give the thermodynamic summary output:

.. literalinclude:: gold.txt


Background
==========

**Ideal gas.** The conversion of electronic structure calculations to
thermodynamic properties in the ideal-gas limit is well documented; see, for
example, Chapter 10 of Cramer, 2004. The key equations used in the
:class:`IdealGasThermo` class are summarized here.

   C.J. Cramer. *Essentials of Computational Chemistry*, Second Edition.
   Wiley, 2004.

The ideal-gas enthalpy is calculated from extrapolation of the energy at 0 K
to the relevant temperature (for an ideal gas, the enthalpy is not a function
of pressure):

.. math ::
   H(T) = E_\text{elec} + E_\text{ZPE} + \int_0^\text{T} C_P \, \text{d}T

where the first two terms are the electronic energy and the zero-point energy,
and the integral is over the constant-pressure heat capacity. The heat
capacity is separable into translational, rotational, vibrational, and
electronic parts (plus a term of :math:`k_\text{B}` to switch from
constant-volume to constant-pressure):

.. math ::
   C_P = k_\text{B} + C_{V\text{,trans}} + C_{V\text{,rot}} + C_{V\text{,vib}} + C_{V\text{,elec}}

The translational heat capacity is 3/2 :math:`k_\text{B}` for a 3-dimensional
gas. The rotational heat capacity is 0 for a monatomic species,
:math:`k_\text{B}` for a linear molecule, and 3/2 :math:`k_\text{B}` for a
nonlinear molecule. In this module, the electronic component of the heat
capacity is assumed to be 0. The vibrational heat capacity contains
:math:`3N-6` degrees of freedom for nonlinear molecules and :math:`3N-5`
degrees of freedom for linear molecules (where :math:`N` is the number of
atoms). The integrated form of the vibrational heat capacity is:

.. math ::
   \int_0^T C_{V,\text{vib}} \text{d}T = \sum_i^\text{vib DOF}
   \frac{\epsilon_i}{e^{\epsilon_i / k_\text{B} T} - 1 }

where :math:`\epsilon_i` are the energies associated with the vibrational
frequencies, :math:`\epsilon_i = h \omega_i`.

The ideal gas entropy can be calculated as a function of temperature and
pressure as:

.. math ::
   S(T,P) &= S(T,P^\circ) - k_\text{B} \ln \frac{P}{P^\circ} \\
          &= S_\text{trans} + S_\text{rot} + S_\text{elec} + S_\text{vib} - k_\text{B} \ln \frac{P}{P^\circ}

where the translational, rotational, electronic, and vibrational components
are calculated as below. (Note that the translational component also includes
components from the Stirling approximation, and that the vibrational degrees
of freedom are enumerated the same as in the above.)

.. math ::
   S_\text{trans} = k_\text{B} \left\{ \ln \left[ \left(
   \frac{2 \pi M k_\text{B} T}{h^2} \right)^{3/2}
   \frac{k_\text{B} T}{P^\circ} \right] + \frac{5}{2} \right\}

.. math ::
   S_\text{rot} = \left\{  \begin{array}{ll}
   0 & \text{, if monatomic} \\
   k_\text{B} \left[ \ln \left( \frac{8\pi^2 I k_\text{B}T}{\sigma h^2}\right) + 1 \right] & \text{, if linear} \\
   k_\text{B} \left\{ \ln \left[ \frac{\sqrt{\pi I_\text{A} I_\text{B} I_\text{C}}}{\sigma} \left(\frac{8\pi^2 k_\text{B} T}{h^2}\right)^{3/2}\right] + \frac{3}{2} \right\} & \text{, if nonlinear} \\
   \end{array}
   \right.

.. math ::
   S_\text{vib} = k_\text{B} \sum_i^\text{vib DOF}
   \left[ \frac{\epsilon_i}{k_\text{B}T\left(e^{\epsilon_i/k_\text{B}T}-1\right)} - \ln \left( 1 - e^{-\epsilon_i/k_\text{B}T} \right)\right]

.. math ::
   S_\text{elec} = k_\text{B} \ln \left[
   2 \times \left(\text{spin multiplicity}\right) + 1\right]

:math:`I_\text{A}` through :math:`I_\text{C}` are the three principle moments
of inertia for a non-linear molecule. :math:`I` is the degenerate moment of
inertia for a linear molecule. :math:`\sigma` is the symmetry number of the
molecule.

The ideal-gas Gibbs free energy is then just calculated from the combination
of the enthalpy and entropy:

.. math ::
   G(T,P) = H(T) - T\, S(T,P)

**Hindered translator / hindered rotor.** The conversion of electronic
structure calculations to thermodynamic properties in the hindered
translator / hindered rotor model was developed for adsorbates on close packed
surfaces and is documented by Sprowl, Campbell, and Arnadottir, 2016. The key
equations used in the :class:`HinderedThermo` class are summarized here.

   L.H. Sprowl, C.T. Campbell, and L. Arnadottir. Hindered Translator and
   Hindered Rotor Models for Adsorbates: Partition Functions and Entropies.
   *J. Phys. Chem. C*, **2016**, 120 (18), pp 9719-9731.

   L.H. Sprowl, C.T. Campbell, and L. Arnadottir. Correction to "Hindered
   Translator and Hindered Rotor Models for Adsorbates: Partition Functions and
   Entropies". *J. Phys. Chem. C*, **2017**, 121 (17), pp 9655-9655.

   C.T. Campbell, L.H. Sprowl, and L. Arnadottir. Equilibrium Constants and
   Rate Constants for Adsorbates: Two-Dimensional (2D) Ideal Gas, 2D Ideal
   Lattice Gas, and Ideal Hindered Translator Models. *J. Phys. Chem. C*,
   **2016**, 120 (19), pp 10283-10297.

The :math:`3N-3` largest vibrational frequencies are used to calculate the
vibrational contributions to the internal energy and the entropy. The
remaining three degrees of freedom are calculated from two translational
contributions and one rotational contribution of the adsorbate. The energy
barriers for the adsorbate to translate and rotate on a close packed surface
are used to calculate the translational and rotational frequencies,
respectively. From the translational and rotational frequencies, the
translational and rotational contributions to the internal energy and the
entropy of the adsorbate are determined. The calculation of the translational
frequency is:

.. math ::
   \nu_{trans} = \sqrt{\frac{W_{trans}}{2mA}}

where :math:`W_{trans}` is the translational energy barrier, :math:`m` is the
mass of the adsorbate, and :math:`A` is the area per surface atom, or the
inverse of the surface site density. The rotational frequency is calculated
as:

.. math ::
   \nu_{rot} = \frac{1}{2\pi}\sqrt{\frac{n^2W_{rot}}{2I}}

where :math:`W_{rot}` is the rotational energy barrier, :math:`n` is the
number of equivalent energy minima in a full rotation of the adsorbate, and
:math:`I` is the reduced moment of inertia of the adsorbate about its surface
bond. Two variables are now introduced, a unitless temperature

.. math ::
   T_i = \frac{kT}{h\nu_i}

and a unitless energy barrier

.. math ::
   r_i = \frac{W_i}{h\nu_i}

to ease the internal energy and entropy calculations.

The internal energy of the adsorbate is calculated as:

.. math ::
   U(T) = E_\text{elec} + E_\text{ZPE} + E_\text{trans} + E_\text{rot} + E_\text{vib}

where :math:`E_{trans}` and :math:`E_{rot}` are:

.. math ::
   E_i = k_\text{B}T \left( \frac{1/T_i}{\exp\left[1/T_i\right]-1} -\frac{1}{2} - \frac{1}{\left(2+16r_i\right)T_i} + \frac{r_i}{2T_i} \left( 1 - \frac{\text{I}_1\left[r_i/2T_i\right]}{\text{I}_0\left[r_i/2T_i\right]}\right) \right)

where :math:`I_{n}` is the nth-order modified Bessel function of the first
kind. Similarly for the harmonic limit, :math:`E_{vib}` is:

.. math ::
   E_\text{vib} = k_\text{B}T \sum_i^\text{3N-3} \left( \frac{1/T_i}{\exp\left[1/T_i\right]-1} \right)

The entropy of the adsorbate is calculated as:

.. math ::
   S = S_\text{trans} + S_\text{rot} + S_\text{vib} + S_\text{con}

where :math:`S_{trans}` and :math:`S_{rot}` are:

.. math ::
   S_i = k_\text{B} \left( \frac{1/T_i}{\exp\left[1/T_i\right]-1} - \ln \left[ 1 - \exp\left[-\frac{1}{T_i}\right]\right] - \frac{1}{2} - \frac{r_i}{2T_i}\frac{\text{I}_1\left[r_i/2T_i\right]}{\text{I}_0\left[r_i/2T_i\right]} + \ln\left[\left(\frac{\pi r_i}{T_i}\right)^{1/2}\text{I}_0\left[\frac{r_i}{2T_i}\right]\right] \right)

and :math:`S_{vib}` is:

.. math ::
   S_\text{vib} = k_\text{B} \sum_i^\text{3N-3} \left( \frac{1/T_i}{\exp\left[1/T_i\right]-1} - \ln \left[ 1 - \exp\left[-\frac{1}{T_i}\right]\right] \right)

:math:`S_{con}` is a concentration related entropy and is calculated as:

.. math ::
   S_\text{con} = k_\text{B} \left( 1 - \ln\left[A\left(\frac{N}{A}\right)^0\right] \right)

where

.. math ::
   \left(\frac{N}{A}\right)^0 = e^{1/3}\left(\frac{N_A \text{ 1 bar}}{RT}\right)

The Helmholtz free energy is calculated as:

.. math ::
   F(T) = U(T) - T\, S(T)

If the user assumes that the :math:`pV` term in :math:`H = U + pV` is
negligible, then the Helmholtz free energy can be used to approximate the
Gibbs free energy, as :math:`G = F + pV`.

**Harmonic limit.** The conversion of electronic structure calculation
information into thermodynamic properties is less established for adsorbates.
However, the simplest approach often taken is to treat all :math:`3N` degrees
of freedom of the adsorbate harmonically since the adsorbate often has no
real translational or rotational degrees of freedom. This is the approach
implemented in the :class:`HarmonicThermo` class. Thus,
the internal energy and entropy of the adsorbate are calculated as

.. math ::
   U(T) = E_\text{elec} + E_\text{ZPE} + \sum_i^\text{harm DOF} \frac{\epsilon_i}{e^{\epsilon_i / k_\text{B} T} - 1 }

.. math ::
   S = k_\text{B} \sum_i^\text{harm DOF}
   \left[ \frac{\epsilon_i}{k_\text{B}T\left(e^{\epsilon_i/k_\text{B}T}-1\right)} - \ln \left( 1 - e^{-\epsilon_i/k_\text{B}T} \right)\right]

and the Helmholtz free energy is calculated as

.. math ::
   F(T) = U(T) - T\, S(T)

In this case, the number of harmonic energies (:math:`\epsilon_i`) used in
the summation is generally :math:`3N`, where :math:`N` is the number of atoms
in the adsorbate. If the user assumes that the :math:`pV` term in
:math:`H = U + pV` is negligible, then the Helmholtz free energy can be used
to approximate the Gibbs free energy, as :math:`G = F + pV`.

**Crystalline solid** The derivation of the partition function for a
crystalline solid is fairly straight-forward and can be found, for example,
in Chapter 11 of McQuarrie, 2000.

   D.A. McQuarrie. *Statistical Mechanics*. University Science Books, 2000.

The treatment implemented in the :class:`CrystalThermo` class depends on
introducing normal coordinates to the entire crystal and treating each atom
in the lattice as an independent harmonic oscillator. This yields the
partition function

.. math ::
   Z = \prod_{j=1}^\text{3N} \left( \frac{e^{-\frac{1}{2}\epsilon_j/k_\text{B}T}}{1 - e^{-\epsilon_j/k_\text{B}T}} \right) e^{-E_\text{elec} / k_\mathrm{B}T}

where :math:`\epsilon_j` are the :math:`3N` vibrational energy levels and
:math:`E_\text{elec}` is the electronic energy of the crystalline solid.
Now, taking the logarithm of the partition function and replacing the
resulting sum with an integral (assuming that the energy level spacing
is essentially continuous) gives

.. math ::
   -\ln Z = E_\text{elec}/k_\text{B}T + \int_0^\infty \left[ \ln \left( 1 - e^{-\epsilon/k_\text{B}T} \right) + \frac{\epsilon}{2 k_\text{B} T} \right]\sigma (\epsilon) \text{d}\epsilon

Here :math:`\sigma (\epsilon)` represents the degeneracy or phonon density of
states as a function of vibrational energy. Once this function has been
determined (i.e. using the :mod:`ase.phonons` module), it is a simple matter
to calculate the canonical ensemble thermodynamic quantities; namely the
internal energy, the entropy and the Helmholtz free energy.

.. math ::
   U(T) &= -\left( \frac{\partial \ln Z}{\partial \frac{1}{k_\text{B}T} } \right)_\text{N,V} \\
        &= E_\text{elec} + \int_0^\infty \left[ \frac{\epsilon}{e^{\epsilon/k_\text{B}T} - 1} + \frac{\epsilon}{2} \right]\sigma (\epsilon) \text{d}\epsilon

.. math ::
   S(T) &= \frac{U}{T} + k_\text{B} \ln Z \\
        &= \int_0^\infty \left[ \frac{\epsilon}{T} \frac{1}{e^{\epsilon/k_\text{B}T} - 1} - k_\text{B} \ln \left(1 - e^{-\epsilon/k_\text{B}T} \right) \right]\sigma (\epsilon) \text{d}\epsilon

.. math ::
   F(T) = U(T) - T\, S(T,P)
