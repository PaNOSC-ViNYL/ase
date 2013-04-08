from ase.cli.command import Command
from ase.io import write


class WriteCommand(Command):
    def add_parser(self, subparser):
        parser = subparser.add_parser('write', help='write ...')
        parser.add_argument('filename', nargs='?')
        
    def run(self, atoms, name):
        print self.args
        filename = self.args.filename or '.traj'
        if 0:#filename[0] == '.':
            filename = name + tag+filename
        #write(filename, atoms)
        write(name + '.json', atoms)
