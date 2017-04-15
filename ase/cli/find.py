import multiprocessing as mp
import os
import os.path as op

from ase.io import read
from ase.io.formats import filetype
# from ase.db.core import parse_query, check_xxx


class CLICommand:
    short_description = 'Find files with atoms in.'

    @staticmethod
    def add_arguments(parser):
        parser.add_argument('folder')
        parser.add_argument('query')
        parser.add_argument('-v', '--verbose', action='store_true')

    @staticmethod
    def run(args):
        main(args)


def main(args):
    # keys, cmps = parse_query(...)

    paths = mp.Queue()
    results = mp.Queue()

    for dirpath, dirnames, filenames in os.walk(args.folder):
        for name in filenames:
            path = op.join(dirpath, name)
            paths.put(path)
        # Skip .git and friends:
        dirnames[:] = (name for name in dirnames if name[0] != '.')

    with mp.Pool() as pool:
        pool.apply(check, [paths, results])

    print(results)


def check(paths, results):
    for path in iter(paths.get, None):
        try:
            filetype(path)
        except ValueError:
            continue
        # atoms = read(path)
        if 1:  # check(atoms):
            results.put(path)
