from ase.utils import prnt


class Command:
    def __init__(self, logfile, args, runfunc=None):
        self.logfile = logfile
        self.args = args
        self.runfunc = runfunc

    def log(self, *args, **kwargs):
        prnt(file=self.logfile, *args, **kwargs)

    @classmethod
    def add_parser(cls, subparser):
        pass

    def run(self, atoms, name):
        self.runfunc(atoms, name, self.args)
