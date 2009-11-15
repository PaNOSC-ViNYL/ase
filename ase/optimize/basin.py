import numpy as np

from ase.optimize import Dynamics
from ase.optimize.fire import FIRE
from ase.units import kB
from ase.parallel import world

class BasinHopping(Dynamics):
    """Basin hopping algorythm.

    After Wales and Doye, J. Phys. Chem. A, vol 101 (1997) 5111-5116"""

    def __init__(self, atoms,
                 temperature=100 * kB,
                 optimizer=FIRE,
                 fmax=0.1,
                 dr=.1,
                 logfile='-', 
                 trajectory=None,
                 adjust_cm=True):
        Dynamics.__init__(self, atoms, logfile, trajectory)
        self.kT = temperature
        self.optimizer = optimizer
        self.fmax = fmax
        self.dr = dr
        if adjust_cm:
            self.cm = atoms.get_center_of_mass()
        else:
            self.cm = None
        
    def run(self, steps):
        """Hop the basins for defined number of steps."""

        # initialize
        atoms = self.atoms
        Eo = self.energy(atoms)
        ro = atoms.get_positions()
        self.Emin = Eo
        self.rmin = ro

        for step in range(steps):
            rn = self.move(ro)
 
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

    def move(self, ro):
        atoms = self.atoms
        # displace coordinates
        disp = np.random.uniform(-1., 1., (len(atoms), 3))
        rn = ro + self.dr * disp
        atoms.set_positions(rn)
        if self.cm is not None:
            cm = atoms.get_center_of_mass()
            atoms.translate(self.cm - cm)
        rn = atoms.get_positions()
        if world is not None:
            world.broadcast(rn, 0)
        atoms.set_positions(rn)
        return atoms.get_positions()

    def get_minimum(self):
        self.atoms.set_positions(self.rmin)
        return self.Emin, self.atoms

    def energy(self, atoms):
        """Return the energy of the nearest local minimum."""
        opt = self.optimizer(atoms)
        opt.run(fmax=self.fmax)
        return atoms.get_potential_energy()
        
