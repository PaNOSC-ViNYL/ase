from ase.utils import prnt
from ase.tasks.io import read_json


class Command:
    def __init__(self, logfile, args):
        self.logfile = logfile
        self.args = args

    def log(self, *args, **kwargs):
        prnt(file=self.logfile, *args, **kwargs)

    @classmethod
    def add_parser(cls, subparser):
        pass

    def run(self, atoms, name):
        pass

    def finalize(self):
        pass
    
    def read(self):
        filename = 'asec' + self.args.tag + '.json'
        return read_json(filename)
