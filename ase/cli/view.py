from ase.asec.command import Command
from ase.visualize import view


class ViewCommand(Command):
    def add_parser(self, subparser):
        parser = subparser.add_parser('view', help='View atoms with ase-gui')
        
    def run(self, atoms, name):
        view(atoms)
