from ase.asec.run import RunCommand
from ase.constraints import FixAtoms
from ase.optimize import LBFGS


class OptimizeCommand(RunCommand):
    @classmethod
    def add_parser(cls, subparser):
        parser = subparser.add_parser('optimize', help='relax ...')
        cls.add_arguments(parser)
        
    @classmethod
    def add_arguments(cls, parser):
        RunCommand.add_arguments(parser)
        parser.add_argument('-f', '--maximum-force', default=0.05, type=float,
                            help='Relax internal coordinates.')
        parser.add_argument('--constrain-tags',
                            metavar='T1,T2,...',
                            help='Constrain atoms with tags T1, T2, ...')
        parser.add_argument('--maximum-stress', type=float,
                            help='Relax unit-cell.')

    def calculate(self, atoms, name):
        args = self.args
        if args.constrain_tags:
            tags = [int(t) for t in args.constrain_tags.split(',')]
            mask = [t in tags for t in atoms.get_tags()]
            atoms.constraints = FixAtoms(mask=mask)
        optimizer = LBFGS(atoms,
                          trajectory=self.get_filename(name, 'traj'),
                          logfile=self.logfile)
        optimizer.run(fmax=args.maximum_force)
        data = RunCommand.calculate(self, atoms, name)

        if hasattr(optimizer, 'force_calls'):
            data['force calls'] = optimizer.force_calls

        return data

