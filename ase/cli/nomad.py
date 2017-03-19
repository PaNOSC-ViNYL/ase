import os
import os.path as op
import subprocess

description = 'Upload files to NOMAD'


def add_arguments(parser):
    parser.add_argument('folders', nargs='+')
    parser.add_argument('-t', '--trajectory')
    parser.add_argument('-k', '--key')
    parser.add_argument('-n', '--do-not-save-key', action='store_true')


def main(args):
    dotase = op.expanduser('~/.ase')
    keyfile = op.join(dotase, 'nomad-key')
    if args.key:
        key = args.key
    else:
        try:
            with open(keyfile) as fd:
                key = fd.readline()
        except FileNotFoundError:
            print()
            return

    subprocess.check_call('echo {} {} | cat'.format(' '.join(args.folders),
                                                    key),
                          shell=True)
    if args.key and not args.do_not_save_key:
        msg = 'Should I save your key for later use? (yes/No): '
        if input(msg).lower() != 'yes':
            if not op.isdir(dotase):
                os.mkdir(dotase)
            with open(keyfile, 'w') as fd:
                print(key, file=fd)
            os.chmod(keyfile, 0o600)
