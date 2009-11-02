import numpy as np

from ase.optimize import Dynamics
from ase.optimize.fire import FIRE
from ase.units import kB

class BasinHopping(Dynamics):
    """Basin hopping algorythm.

    After Wales and Doye, J. Phys. Chem. A, vol 101 (1997) 5111-5116"""

    def __init__(self, atoms,
                 temperature=100 * kB,
                 optimizer=FIRE,
                 fmax=0.1,
                 dr=.1,
                 logfile='-', trajectory=None):
        Dynamics.__init__(self, atoms, logfile, trajectory)
        self.kT = temperature
        self.optimizer = optimizer
        self.fmax = fmax
        self.dr = dr
        
    def run(self, steps):
        """Hop the basins for defined number of steps."""

        # initialize
        atoms = self.atoms
        Eo = self.energy(atoms)
        ro = atoms.get_positions()
        self.Emin = Eo
        self.rmin = ro

        for step in range(steps):
            # displace coordinates
            disp = np.random.uniform(-1., 1., (len(atoms), 3))
            rn = ro + self.dr * disp
            atoms.set_positions(rn)
 
            En = self.energy(atoms)
            if En < self.Emin:
                # new minimum found
                self.Emin = En
                self.rmin = atoms.get_positions()
                self.call_observers()

            accept = np.exp((Eo - En) / self.kT) > np.random.uniform()
            if accept:
                ro = rn
                Eo = En

            if self.logfile is not None:
                name = self.__class__.__name__
                self.logfile.write('%s: step %d, energy %15.6f, emin %15.6f\n'
                                   % (name, step, En, self.Emin))
                self.logfile.flush()

    def get_minimum(self):
        self.atoms.set_positions(self.rmin)
        return self.Emin, self.atoms

    def energy(self, atoms):
        """Return the energy of the nearest local minimum."""
        opt = self.optimizer(atoms, logfile=None)
        opt.run(fmax=self.fmax)
        return atoms.get_potential_energy()
        
