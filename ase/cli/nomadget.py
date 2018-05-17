from __future__ import print_function
import json


class CLICommand:
    short_description = 'Get calculation from NOMAD and write json to stdout'

    @staticmethod
    def add_arguments(p):
        p.add_argument('uri', metavar='nmd://<hash>',
                       help='URI to get')

    @staticmethod
    def run(args):
        from ase.nomad import download
        calculation = download(args.uri)
        print(json.dumps(calculation.dct))
