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
    query = 42

    N = mp.cpu_count()
    paths = mp.Queue()
    results = mp.Queue()

    pool = mp.Pool(N, initialize_target, [paths, results])
    result = pool.map_async(target, [query] * N)

    for path in allpaths(args):
        paths.put(path)
    for _ in range(N):
        paths.put(None)

    result.get()
    pool.terminate()

    while not results.empty():
        print(results.get())


def allpaths(args):
    for dirpath, dirnames, filenames in os.walk(args.folder):
        for name in filenames:
            if name.endswith('.py'):
                continue
            if name.endswith('py'):
                continue
            path = op.join(dirpath, name)
            yield path
        # Skip .git, __pycache__ and friends:
        dirnames[:] = (name for name in dirnames if name[0] not in '._')


def target(query):
    for path in iter(target.paths.get, None):
        if check(path, query):
            target.results.put(path)


def check(path, query):
    format = filetype(path, guess=False)
    if format is None:
        return False
    #atoms = read(path, format=format)
    return True


def initialize_target(paths, results):
    target.paths = paths
    target.results = results
