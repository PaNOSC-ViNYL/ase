from ase.asec.run import RunCommand
from ase.constraints import FixAtoms


class OptimizeCommand(RunCommand):
    @classmethod
    def add_parser(cls, subparser):
        parser = subparser.add_parser('relax', help='relax ...')
        cls.add_arguments(parser)
        
    @classmethod
    def add_arguments(cls, parser):
        RunCommand.add_arguments(parser)
        parser.add_argument('-R', '--relax', metavar='FMAX[,OPTIMIZER]',
                            help='Relax internal coordinates using '
                            'OPTIMIZER algorithm. The OPTIMIZER keyword is '
                            'optional, and if omitted LBFGS is used by default.')
        parser.add_argument('--relaxsteps', type=int,
                            metavar='steps',
                            help='Limit the number of optimizer steps.')
        parser.add_argument('--constrain-tags',
                            metavar='T1,T2,...',
                            help='Constrain atoms with tags T1, T2, ...')
        parser.add_argument('--srelax', metavar='SFMAX[,SOPTIMIZER]',
                        help='Relax cell by minimizing stress using StranFilter '
                        'with SOPTIMIZER algorithm. The SOPTIMIZER keyword is '
                        'optional, and if omitted BFGS is used by default.')

    def run(self, atoms, name):

        if opts.relax:
            if len(opts.relax.split(',')) > 1:
                self.fmax, self.optimizer = opts.relax.split(',')
            else:
                self.fmax = opts.relax
                self.optimizer = 'LBFGS'
            self.fmax = float(self.fmax)

        if opts.relaxsteps is not None:
            self.steps = int(opts.relaxsteps)
        else:
            # yes, the default number of ASE optimizer steps
            # ase/optimize/optimize.py
            self.steps = 100000000

        if opts.constrain_tags:
            self.constrain_tags = [int(t)
                                   for t in opts.constrain_tags.split(',')]

        mask = [t in self.constrain_tags for t in atoms.get_tags()]
        if mask:
            constrain = FixAtoms(mask=mask)
            atoms.constraints = [constrain]

        optstr = "ase.optimize." + self.optimizer
        optstr += "(atoms, "
        if trajectory is not None:
            # allow existing trajectory to append new optimizations
            # or create new one if trajectory is str
            optstr += "trajectory=trajectory,"
        else:
            # create new trajectory with default name
            optstr += "trajectory=self.get_filename(name, 'traj'),"
        optstr += "logfile=self.logfile)"
        optimizer = eval(optstr)
        try:
            # handle scipy optimizers who raise Converged when done
            from ase.optimize import Converged
            try:
                optimizer.run(fmax=self.fmax, steps=self.steps)
            except Converged:
                pass
        except ImportError:
            optimizer.run(fmax=self.fmax, steps=self.steps)
        # StrainFilter optimizer steps
        steps = optimizer.get_number_of_steps()
        if data.get('optimizer steps', None) is None:
            data['optimizer steps'] = steps
        else:
            data['optimizer steps'] += steps
        # optimizer force calls
        if hasattr(optimizer, 'force_calls'):
            calls = optimizer.force_calls
        else:
            calls = steps
        if data.get('optimizer force calls', None) is None:
            data['optimizer force calls'] = calls
        else:
            data['optimizer force calls'] += calls

        if self.fmax is not None:
            self.optimize(name, atoms, data)

            e = atoms.get_potential_energy()

            data['relaxed energy'] = e

        return data

