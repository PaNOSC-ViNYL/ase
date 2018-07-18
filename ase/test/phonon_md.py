import numpy as np
from numpy.random import RandomState
import ase.units as u
from ase.md.verlet import VelocityVerlet
from ase.phonons import Phonons
from ase.data import atomic_numbers
from ase.optimize import FIRE
#from asap3 import EMT
from ase.calculators.emt import EMT
from ase.build import bulk
from ase.md.velocitydistribution import PhononHarmonics

rng = RandomState(17)

atoms = bulk('Pd')
atoms *= (3, 3, 3)
avail = [atomic_numbers[sym]
         for sym in ['Ni', 'Cu', 'Pd', 'Ag', 'Pt', 'Au']]
atoms.numbers[:] = rng.choice(avail, size=len(atoms))
atoms.calc = EMT()

opt = FIRE(atoms, trajectory='relax.traj')
opt.run(fmax=0.001)
positions0 = atoms.positions.copy()

phonons = Phonons(atoms, EMT(), supercell=(1, 1, 1), delta=0.05)

try:
    phonons.run()
    phonons.read()  # Why all this boilerplate?
finally:
    phonons.clean()
matrices = phonons.get_force_constant()

K = matrices[0]
T = 100 * u.kB

atoms.calc = EMT()
Epotref = atoms.get_potential_energy()

temps = []
Epots = []
Ekins = []
Etots = []


for i in range(100):
    PhononHarmonics(atoms, T, K, rng=np.random.RandomState(888 + i))

    Epot = atoms.get_potential_energy() - Epotref
    Ekin = atoms.get_kinetic_energy()
    Ekins.append(Ekin)
    Epots.append(Epot)
    Etots.append(Ekin + Epot)
    temps.append(atoms.get_temperature())

    atoms.positions[:] = positions0

    # The commented code would produce displacements/velocities
    # resolved over phonon modes if we borrow some expressions
    # from the function.  Each mode should contribute on average
    # equally to both Epot and Ekin/temperature
    #
    # atoms1.calc = EMT()
    # atoms1 = atoms.copy()
    # v_ac = np.zeros_like(positions0)
    # D_acs, V_acs = ...
    # for s in range(V_acs.shape[2]):
    #     atoms1.positions += D_acs[:, :, s]
    #     v_ac += V_acs[:, :, s]
    #     atoms1.set_velocities(v_ac)
    #     X1.append(atoms1.get_potential_energy() - Epotref)
    #     X2.append(atoms1.get_kinetic_energy())

    print('energies', Epot, Ekin, Epot + Ekin)



Epotmean = np.mean(Epots)
Ekinmean = np.mean(Ekins)
Tmean = np.mean(temps)
Terr = abs(Tmean - T / u.kB)
relative_imbalance = abs(Epotmean - Ekinmean) / (Epotmean + Ekinmean)


assert Epotmean - Ekinmean

print('epotmean', Epotmean)
print('ekinmean', Ekinmean)
print('rel imbalance', relative_imbalance)
print('Tmean', Tmean, 'Tref', T / u.kB, 'err', Terr)

assert Terr < 4.0  # error in Kelvin for instantaneous velocity
assert relative_imbalance < 0.02  # Epot == Ekin give or take 2 %

if 0:
    import matplotlib.pyplot as plt
    I = np.arange(len(Epots))
    plt.plot(I, Epots, 'o', label='pot')
    plt.plot(I, Ekins, 'o', label='kin')
    plt.plot(I, Etots, 'o', label='tot')
    plt.show()

