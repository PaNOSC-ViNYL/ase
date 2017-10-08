from ase.io import read


class CLICommand:
    short_description = 'Convert file formats'

    @staticmethod
    def add_arguments(parser):
        parser.add_argument('-v', '--verbose', action='store_true')
        parser.add_argument('input', metavar='input-file')
        parser.add_argument('-i', '--input-format', type=str)
        parser.add_argument('output', metavar='output-file')
        parser.add_argument('-o', '--output-format', type=str)

    @staticmethod
    def run(args):
        if args.verbose:
            print(args.input + " -> " + args.output)

        atoms = read(args.input, format=args.input_format)
        atoms.write(args.output, format=args.output_format)
