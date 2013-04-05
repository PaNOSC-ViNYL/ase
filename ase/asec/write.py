from ase.asec.command import Command
from ase.io import write


class WriteCommand(Command):
    def add_parser(self, subparser):
        parser = subparser.add_parser('write', help='write ...')
        parser.add_argument('filename', nargs='?')
        
    def run(self, atoms, name):
        filename = self.args.filename or '.traj'
        if filename[0] == '.':
            filename = name + self.args.tag + filename
        write(filename, atoms)
