import numpy as np
from numpy.random import RandomState
import ase.units as u
from ase.md.verlet import VelocityVerlet
from ase.phonons import Phonons
from ase.data import atomic_numbers
from ase.optimize import FIRE
from asap3 import EMT
#from ase.calculators.emt import EMT
from ase.build import bulk
from ase.md.velocitydistribution import PhononHarmonics

rng = RandomState(17)

atoms = bulk('Pd')
atoms *= (4, 4, 4)
avail = [atomic_numbers[sym]
         for sym in ['Ni', 'Cu', 'Pd', 'Ag', 'Pt', 'Au']]
atoms.numbers[:] = rng.choice(avail, size=len(atoms))
atoms.calc = EMT()

opt = FIRE(atoms, trajectory='relax.traj')
opt.run(fmax=0.001)

#atoms.set_masses()
#print(atoms.get_masses())

phonons = Phonons(atoms, EMT(), supercell=(1, 1, 1), delta=0.01)

try:
    phonons.run()  # ugly
    phonons.read(acoustic=True)  # wtf!
finally:
    phonons.clean()
matrices = phonons.get_force_constant()
K = matrices[0]

T = 100 * u.kB

atoms0 = atoms.copy()
atoms0.calc = EMT()
Eref = atoms0.get_potential_energy()

avgs = []
Epots = []
Ekins = []


for i in range(50):
    atoms = atoms0.copy()
    atoms.calc = EMT()
    PhononHarmonics(atoms, T, K, rng=np.random.RandomState(888 + i))
    #dyn = Langevin(atoms, 4 * u.fs, T, 0.01, logfile='-')
    Epot = atoms.get_potential_energy()
    Ekin = atoms.get_kinetic_energy()
    Ekins.append(Ekin)
    Epots.append(Epot - Eref)
    #print(Epot, Ekin)

    #dyn = VelocityVerlet(atoms, 8 * u.fs)
    #temps = []
    #for i in range(100):
    #    f = atoms.get_forces()
    #    tnow = atoms.get_temperature()
        #if i % 100 == 0:
    #    temps.append(tnow)
    #    dyn.step(f)

    #avg = sum(temps) / len(temps)
    #print('average temp', avg)
    #avgs.append(avg)

#print('average average', np.mean(avgs))
    
#dyn.run(1000)
mean1 = np.mean(Epots)
mean2 = np.mean(Ekins)
print(mean1)
print(mean2)
#print(np.
