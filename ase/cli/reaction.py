import numpy as np

from ase.cli.command import Command


class ReactionCommand(Command):
    def add_parser(self, subparser):
        parser = subparser.add_parser('reaction', help='reaction ...')

    def __init__(self):
        self.count = {}

    def run(self, atoms, name):
        from collections import Counter
        self.count[name] = Counter(atoms.numbers)
    
    def finalize(self):
        numbers = set()
        for c in self.count.values():
            numbers.update(c.keys())
        print numbers
        a = []
        for name in self.args.names:
            a.append([self.count[name].get(Z, 0) for Z in numbers])
            
        print a
        u, s, v = np.linalg.svd(a)
        coefs = u[:, -1] / u[0, -1]
        energy = 0.0
        for c, name in zip(coefs, self.args.names):
            e = self.data[name]['energy']
            energy += c * e
            print '%10.5f %f %s' % (c, e, name)
        print energy, 'ev'
