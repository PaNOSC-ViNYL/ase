from ase.io import read, write


class CLICommand:
    short_description = 'Convert between file formats'
    description = 'Convert between file formats.  Use "-" for stdin/stdout.'

    @staticmethod
    def add_arguments(parser):
        add = parser.add_argument
        add('-v', '--verbose', action='store_true')
        add('input', nargs='+', metavar='input-file')
        add('-i', '--input-format')
        add('output', metavar='output-file')
        add('-o', '--output-format')
        add('-n', '--image-number',
            default=':', metavar='NUMBER',
            help='Pick image(s) from trajectory.  NUMBER can be a '
            'single number (use a negative number to count from '
            'the back) or a range: start:stop:step, where the '
            '":step" part can be left out - default values are '
            '0:nimages:1.')

    @staticmethod
    def run(args):
        if args.verbose:
            print(', '.join(args.input), '->', args.output)

        configs = []
        for filename in args.input:
            atoms = read(filename, args.image_number, format=args.input_format)
            if isinstance(atoms, list):
                configs.extend(atoms)
            else:
                configs.append(atoms)

        write(args.output, configs, format=args.output_format)
