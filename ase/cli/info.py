import platform
import sys

from ase.utils import import_module
from ase.io.formats import filetype
from ase.io.ulm import print_ulm_info
from ase.io.pickletrajectory import print_trajectory_info
from ase.io.bundletrajectory import print_bundletrajectory_info

description = 'Print summary of information from trajectory files.'


def add_arguments(parser):
    parser.description = description
    parser.add_argument('filenames', nargs='*')


def main(args):
    for f in args.filenames:
        ft = filetype(f)
        print("File type of '{}' appears to be of type '{}'".format(f, ft))
        if ft == 'traj':
            print_ulm_info(f)
        elif ft == 'trj':
            print_trajectory_info(f)
        elif ft == 'bundletrajectory':
            print_bundletrajectory_info(f)
    if not args.filenames:
        print_info()


def print_info():
    versions = [('platform', platform.platform()),
                ('python-' + sys.version.split()[0], sys.executable)]
    for name in ['ase', 'numpy', 'scipy']:
        try:
            module = import_module(name)
        except ImportError:
            versions.append((name, 'no'))
        else:
            versions.append((name + '-' + module.__version__,
                            module.__file__.rsplit('/', 1)[0] + '/'))

    for a, b in versions:
        print('{:16}{}'.format(a, b))
