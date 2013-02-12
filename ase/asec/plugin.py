from ase.asec.command import Command

class PluginCommand(Command):
    def __init__(self, logfile, args, run_function):
        RunCommand.__init__(self, logfile, args)
        self.run_function = run_function

    def run(self, atoms, name, args):
        self.run_function(atoms, name, args)
