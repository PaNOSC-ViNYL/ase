from __future__ import division
import numpy as np

from ase.cli.run import RunCommand
from ase.io.trajectory import PickleTrajectory
from ase.utils.eos import EquationOfState


class EOSCommand(RunCommand):
    def add_parser(self, subparser):
        parser = subparser.add_parser('eos', help='Equation of state ...')
        self.add_arguments(parser)
        
    def add_arguments(self, parser):
        RunCommand.add_arguments(self, parser)
        add = parser.add_argument
        add('--points', default=5, type=int, help='Number of points for fit.')
        add('--eos-type', default='sjeos', help='Selects the type of eos.')

    def calculate(self, atoms, name):
        args = self.args
        
        traj = PickleTrajectory(self.get_filename(name, 'traj'), 'w', atoms)
        eps = 0.01
        strains = np.linspace(1 - eps, 1 + eps, args.points)
        v1 = atoms.get_volume()
        volumes = strains**3 * v1
        energies = []
        cell1 = atoms.cell
        for s in strains:
            atoms.set_cell(cell1 * s, scale_atoms=True)
            energies.append(atoms.get_potential_energy())
            traj.write(atoms)
        traj.close()
        eos = EquationOfState(volumes, energies, args.eos_type)
        v0, e0, B = eos.fit()
        atoms.set_cell(cell1 * (v0 / v1)**(1 / 3), scale_atoms=True)
        data = RunCommand.calculate(self, atoms, name)
        data.update({'volumes': volumes,
                     'energies': energies,
                     'fitted_energy': e0,
                     'fitted_volume': v0,
                     'bulk_modulus': B,
                     'eos_type': args.eos_type})
        return data


EosCommand = EOSCommand
