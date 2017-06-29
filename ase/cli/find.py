from __future__ import print_function
import multiprocessing as mp
import os
import os.path as op
import sys

from ase.io import read
from ase.io.formats import filetype
from ase.db import connect
from ase.db.core import parse_selection
from ase.db.jsondb import JSONDatabase
from ase.db.row import atoms2dict


class CLICommand:
    short_description = 'Find files with atoms in.'

    @staticmethod
    def add_arguments(parser):
        parser.add_argument('folder')
        parser.add_argument('query', nargs='?')
        parser.add_argument('-v', '--verbose', action='store_true')
        parser.add_argument('-i', '--include')
        parser.add_argument('-x', '--exclude')
        parser.add_argument('-j', '--jobs', type=int, default=0)

    @staticmethod
    def run(args):
        main(args)


def main(args):
    query = parse_selection(args.query)
    include = args.include.split(',') if args.include else []
    exclude = args.exclude.split(',') if args.exclude else []
    N = args.jobs or mp.cpu_count()
    if N == 1:
        for path in allpaths(args.folder, include, exclude):
            if check(path, query, args.verbose):
                print(path)
        return

    paths = mp.Queue()

    pool = mp.Pool(N, initialize_target, [paths])

    results = [pool.apply_async(target, (query, args.verbose))
               for _ in range(N)]

    for path in allpaths(args.folder, include, exclude):
        paths.put(path)
    for _ in range(N):
        paths.put('STOP')

    found = []
    for result in results:
        found.extend(result.get())

    pool.terminate()

    for path in sorted(found):
        print(path)


def allpaths(folder, include, exclude):
    exclude += ['.py', '.pyc']
    for dirpath, dirnames, filenames in os.walk(folder):
        for name in filenames:
            if any(name.endswith(ext) for ext in exclude):
                continue
            if include:
                for ext in include:
                    if name.endswith(ext):
                        break
                else:
                    continue
            path = op.join(dirpath, name)
            yield path

        # Skip .git, __pycache__ and friends:
        dirnames[:] = (name for name in dirnames if name[0] not in '._')


def target(query, verbose):
    paths = []
    for path in iter(target.paths.get, 'STOP'):
        if check(path, query, verbose):
            paths.append(path)
    return paths


def check(path, query, verbose):
    try:
        format = filetype(path, guess=False)
    except OSError:
        return False
    if format is None:
        return False
    if format in ['db', 'json']:
        db = connect(path)
    else:
        try:
            atoms = read(path, format=format)
        except Exception as x:
            if verbose:
                print(path + ':', x, file=sys.stderr)
            return False
        db = FakeDB(atoms)

    try:
        for row in db._select(*query):
            return True
    except Exception as x:
        if verbose:
            print(path + ':', x, file=sys.stderr)

    return False


def initialize_target(paths):
    target.paths = paths


class FakeDB(JSONDatabase):
    def __init__(self, atoms):
        self.bigdct = {1: atoms2dict(atoms)}

    def _read_json(self):
        return self.bigdct, [1], 2
