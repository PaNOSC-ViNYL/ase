import optparse

import numpy as np

from ase.lattice import bulk
from ase.tasks.task import OptimizeTask
from ase.data import chemical_symbols
from ase.utils.eos import EquationOfState
from ase.io.trajectory import PickleTrajectory
import ase.units as units


class BulkTask(OptimizeTask):
    taskname = 'bulk'

    def __init__(self, crystal_structure=None, lattice_constant=None,
                 c_over_a=None, cubic=False, orthorhombic=False, fit=False,
                 **kwargs):
        """Bulk task."""

        self.crystal_structure = crystal_structure
        self.lattice_constant = lattice_constant
        self.c_over_a = c_over_a
        self.cubic = cubic
        self.orthorhombic = orthorhombic
        self.fit = fit

        OptimizeTask.__init__(self, **kwargs)
        
        self.summary_header += [('V0', 'Ang^3'),
                                ('B', 'GPa')]

    def build_system(self, name):
        atoms = bulk(name, crystalstructure=self.crystal_structure,
                     a=self.lattice_constant, covera=self.c_over_a,
                     orthorhombic=self.orthorhombic, cubic=self.cubic)

        if self.repeat is not None:
            r = self.repeat.split(',')
            if len(r) == 1:
                r = 3 * r
            atoms = atoms.repeat([int(c) for c in r])

        return atoms
    
    def fit_volume(self, name, atoms):
        cell0 = atoms.get_cell()
        strains = np.linspace(0.98, 1.02, 5)
        energies = []
        traj = PickleTrajectory(self.get_filename(name, '-fit.traj'), 'w')
        for s in strains:
            atoms.set_cell(cell0 * s, scale_atoms=True)
            energies.append(atoms.get_potential_energy())
            traj.write(atoms)

        traj.close()

        data = {'energy': energies[2],
                'strains': strains,
                'energies': energies}

        return data
            
    def calculate(self, name, atoms):
        if self.fit:
            return self.fit_volume(name, atoms)
        else:
            return OptimizeTask.calculate(self, name, atoms)
        
    def analyse(self):
        OptimizeTask.analyse(self)
        for name, data in self.data.items():
            if 'strains' in data:
                atoms = self.create_system(name)
                volumes = data['strains']**3 * atoms.get_volume()
                energies = data['energies']
                eos = EquationOfState(volumes, energies)
                v, e, B = eos.fit()

                self.results[name][1:] = [energies[2] - e, v,
                                          B * 1e24 / units.kJ]
            else:
                self.results[name].extend([None, None])

    def add_options(self, parser):
        OptimizeTask.add_options(self, parser)

        bulk = optparse.OptionGroup(parser, 'Bulk')
        bulk.add_option('-F', '--fit', action='store_true',
                        help='Find optimal volume.')
        bulk.add_option('-x', '--crystal-structure',
                        help='Crystal structure.',
                        choices=['sc', 'fcc', 'bcc', 'diamond', 'hcp',
                                 'rocksalt', 'zincblende'])
        bulk.add_option('-a', '--lattice-constant', type='float',
                        help='Lattice constant in Angstrom.')
        bulk.add_option('--c-over-a', type='float',
                        help='c/a ratio.')
        bulk.add_option('-O', '--orthorhombic', action='store_true',
                        help='Use orthorhombic unit cell.')
        bulk.add_option('-C', '--cubic', action='store_true',
                        help='Use cubic unit cell.')
        bulk.add_option('-r', '--repeat',
                        help='Repeat unit cell.  Use "-r 2" or "-r 2,3,1".')
        parser.add_option_group(bulk)

    def parse(self, opts, args):
        OptimizeTask.parse(self, opts, args)

        self.fit = opts.fit
        self.crystal_structure = opts.crystal_structure
        self.lattice_constant = opts.lattice_constant
        self.c_over_a = opts.c_over_a
        self.orthorhombic = opts.orthorhombic
        self.cubic = opts.cubic
        self.repeat = opts.repeat
