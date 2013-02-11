from ase.asec.run import RunCommand

class PluginCommand(RunCommand):
    @classmethod
    def add_parser(cls, subparser):
        parser = subparser.add_parser('plugin', help='plugin ...')
        cls.add_arguments(parser)

    def run(self):
        if 'run' in namespace:
            self.command = Command(self.logfile, args,
                                   namespace.get('run'))
