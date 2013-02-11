from ase.asec.command import Command


class ReactionCommand(Command):
    @classmethod
    def add_parser(cls, subparser):
        parser = subparser.add_parser('reaction', help='reaction ...')
        cls.add_arguments(parser)

    @classmethod
    def add_arguments(cls, parser):
        parser.add_argument('formula')

    def run(self):
        a=np.array([(1,2,0),(2,4,0),(0,10,10),(0,10,10)])
        u,s,v=np.linalg.svd(a)
        print v[-1]/v[-1,0]
