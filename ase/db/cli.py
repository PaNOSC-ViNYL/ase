import argparse

from ase.db import connect


def main():
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument('name', nargs=1)
    parser.add_argument('query', nargs=1)
    args = parser.parse_args()
    con = connect(args.name[0])
    print args.query
    print list(con.iselect(*args.query))
