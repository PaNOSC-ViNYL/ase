# -*- coding: utf-8 -*-
from __future__ import division, print_function

from ase.optimize.optimize import Optimizer
from ase.units import Ha, Bohr

from berny import Berny as _Berny, Molecule, Logger


class Berny(Optimizer):
    def __init__(self, atoms, restart=None, logfile='-', trajectory=None,
                 master=None, verbosity=-2):
        """Berny optimizer.

        Parameters:

        atoms: Atoms object
            The Atoms object to relax.

        restart: string
            Pickle file used to store internal state. If set, file with
            such a name will be searched and internal state stored will
            be used, if the file exists.

        trajectory: string
            Pickle file used to store trajectory of atomic movement.

        logfile: file object or str
            If *logfile* is a string, a file with that name will be opened.
            Use '-' for stdout.

        master: boolean
            Defaults to None, which causes only rank 0 to save files.  If
            set to true,  this rank will save files.

        verbosity: int
            Controls amount of output, defaults to -2 which does not output
            anything, -1 outputs only energy, 0 outputs all information
        """
        self._restart_data = None  # Optimizer.__init__() may overwrite
        Optimizer.__init__(self, atoms, restart, logfile, trajectory, master)
        geom = Molecule(atoms.get_chemical_symbols(), atoms.positions)
        self._berny = _Berny(
            geom,
            log=Logger(out=self.logfile, verbosity=verbosity),
            debug=True,
            restart=self._restart_data,
            maxsteps=10000000000,  # TODO copied from ase.optimize.Optimizer
            gradientmax=0.,
            gradientrms=0.,
            stepmax=0.,
            steprms=0.,
        )
        next(self._berny)
        # Berny yields the initial geometry the first time because it is
        # typically used as a generator, see berny.optimize()

    def step(self, f):
        energy = self.atoms.get_potential_energy()
        gradients = -self.atoms.get_forces()
        debug = self._berny.send((energy / Ha, gradients / Ha * Bohr))
        self.dump(debug)
        geom = next(self._berny)
        self.atoms.positions[:] = geom.coords

    def read(self):
        self._restart_data = self.load()
