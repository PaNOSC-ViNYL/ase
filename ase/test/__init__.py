from __future__ import print_function
import os
import sys
import subprocess
from multiprocessing import Pool, cpu_count
import tempfile
import unittest
from glob import glob
from distutils.version import LooseVersion
import time
import traceback

import numpy as np

from ase.calculators.calculator import names as calc_names, get_calculator
from ase.parallel import paropen
from ase.utils import devnull
from ase.cli.info import print_info

NotAvailable = unittest.SkipTest

test_calculator_names = []


def require(calcname):
    if calcname not in test_calculator_names:
        raise NotAvailable('use --calculators={0} to enable'.format(calcname))


class CustomTextTestRunner(unittest.TextTestRunner):
    def __init__(self, logname, descriptions=1, verbosity=1):
        self.f = paropen(logname, 'w')
        unittest.TextTestRunner.__init__(self, self.f, descriptions, verbosity)

    def run(self, test):
        stderr_old = sys.stderr
        try:
            sys.stderr = self.f
            testresult = unittest.TextTestRunner.run(self, test)
        finally:
            sys.stderr = stderr_old
        return testresult


class ScriptTestCase(unittest.TestCase):
    def __init__(self, methodname='testfile', filename=None):
        unittest.TestCase.__init__(self, methodname)
        self.filename = filename

    def testfile(self):
        try:
            with open(self.filename) as fd:
                exec(compile(fd.read(), self.filename, 'exec'), {})
        except ImportError as ex:
            module = ex.args[0].split()[-1].replace("'", '').split('.')[0]
            if module in ['scipy', 'matplotlib', 'Scientific', 'lxml',
                          'flask', 'gpaw', 'GPAW', 'netCDF4']:
                raise unittest.SkipTest('no {} module'.format(module))
            else:
                raise

    def id(self):
        return self.filename

    def __str__(self):
        return self.filename.split('test/')[-1]

    def __repr__(self):
        return "ScriptTestCase(filename='%s')" % self.filename


def get_tests(files=None):
    if files:
        files = [os.path.join(__path__[0], f) for f in files]
    else:
        files = glob(__path__[0] + '/*')

    sdirtests = []  # tests from subdirectories: only one level assumed
    tests = []
    for f in files:
        if os.path.isdir(f):
            # add test subdirectories (like calculators)
            sdirtests.extend(glob(f + '/*.py'))
        else:
            # add py files in testdir
            if f.endswith('.py'):
                tests.append(f)
    tests.sort()
    sdirtests.sort()
    tests.extend(sdirtests)  # run test subdirectories at the end
    tests = [test for test in tests if not test.endswith('__.py')]
    return tests


def sleep_forever():
    # KeyboardInterrupts are caught differently by worker threads in
    # Py2/3, and our own except-clause within runtest_subprocess()
    # will not be reached from worker threads that are idle because
    # they have finished executing tasks.  Therefore instead of
    # idling, they must execute this task:
    try:
        while True:
            time.sleep(100000)
    except KeyboardInterrupt:
        pass


def runtest_subprocess(filename):
    # Py2/3 compatibility hack for idle processes:
    if filename == 'sleep_forever':
        sleep_forever()

    sys.stdout = devnull
    basedir = os.path.split(__file__)[0]
    testrelpath = os.path.relpath(filename, basedir)
    data = {'name': testrelpath, 'pid': os.getpid()}
    test = ScriptTestCase(filename=filename)

    # Some tests may write to files with the same name as other tests.
    # Hence, create new subdir for each test:
    cwd = os.getcwd()
    testsubdir = testrelpath.replace('/', '_').replace('.', '_')
    os.mkdir(testsubdir)
    os.chdir(testsubdir)
    t1 = time.time()

    try:
        test.testfile()
    except unittest.SkipTest as ex:
        data['status'] = 'SKIPPED'
        data['error'] = ex
    except AssertionError as ex:
        data['status'] = 'FAIL'
        data['error'] = ex
        data['traceback'] = traceback.format_exc()
    except BaseException as ex:  # Return any error/signal to parent process.
        # Important: In Py2, subprocesses get the Ctrl+C signal whereas in
        # Py3, the parent process does - which is much better.
        # In Py2 it is difficult to kill the processes since the parent
        # process does not know about Ctrl+C.  This is why we
        # return the exception to the parent process and let the
        # parent process smoothly close the pool.
        data['status'] = 'ERROR'
        data['error'] = ex
        data['traceback'] = traceback.format_exc()
    finally:
        t2 = time.time()
        os.chdir(cwd)

    if 'error' not in data:
        data['status'] = 'OK'

    data['time'] = t2 - t1
    return data


def test(verbosity=1, calculators=[],
         stream=sys.stdout, files=None):
    """Main test-runner for ASE."""

    if LooseVersion(np.__version__) >= '1.14':
        # Our doctests need this (spacegroup.py)
        np.set_printoptions(legacy='1.13')

    test_calculator_names.extend(calculators)
    disable_calculators([name for name in calc_names
                         if name not in calculators])

    tests = get_tests(files)

    ts = unittest.TestSuite()

    print_info()

    origcwd = os.getcwd()
    testdir = tempfile.mkdtemp(prefix='ase-test-')
    os.chdir(testdir)
    nprocs = cpu_count()

    print('{:25}{}'.format('test directory', testdir))
    print('{:25}{}'.format('number of processes', nprocs))
    print('{:25}{}'.format('time', time.strftime('%c')))
    print()

    ntests_done = 0
    nproblems = 0
    pool = Pool(nprocs)
    t1 = time.time()

    try:
        for data in pool.imap_unordered(runtest_subprocess, tests
                                        # Hack: See sleep_forever() function.
                                        + ['sleep_forever'] * nprocs):
            print('{pid:5d}: {name:40} {time:2.2f}s {status}'
                  .format(**data))
            if data.get('traceback'):
                print('=' * 78)
                print(data['traceback'].rstrip())
                print('=' * 78)
            ntests_done += 1
            if data['status'] in ['FAIL', 'ERROR']:
                nproblems += 1

            if ntests_done == len(tests):
                # Workers are sleeping:
                break
    except KeyboardInterrupt:
        print()
        print('Interrupted by keyboard')
    finally:
        os.chdir(origcwd)
        pool.terminate()
        pool.join()

        t2 = time.time()
        print('Time elapsed: {:.1f} s'.format(t2 - t1))

    return nproblems


def disable_calculators(names):
    for name in names:
        if name in ['emt', 'lj', 'eam', 'morse', 'tip3p']:
            continue
        try:
            cls = get_calculator(name)
        except ImportError:
            pass
        else:
            def get_mock_init(name):
                def mock_init(obj, *args, **kwargs):
                    raise NotAvailable('use --calculators={0} to enable'
                                       .format(name))
                return mock_init

            def mock_del(obj):
                pass
            cls.__init__ = get_mock_init(name)
            cls.__del__ = mock_del


def cli(command, calculator_name=None):
    if (calculator_name is not None and
        calculator_name not in test_calculator_names):
        return
    proc = subprocess.Popen(' '.join(command.split('\n')),
                            shell=True,
                            stdout=subprocess.PIPE)
    print(proc.stdout.read().decode())
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError('Failed running a shell command.  '
                           'Please set you $PATH environment variable!')


class must_raise:
    """Context manager for checking raising of exceptions."""
    def __init__(self, exception):
        self.exception = exception

    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_value, tb):
        if exc_type is None:
            raise RuntimeError('Failed to fail: ' + str(self.exception))
        return issubclass(exc_type, self.exception)


class CLICommand:
    short_description = 'Test ASE'

    @staticmethod
    def add_arguments(parser):
        parser.add_argument(
            '-c', '--calculators',
            help='Comma-separated list of calculators to test.')
        parser.add_argument('-v', '--verbose', help='verbose mode',
                            action='store_true')
        parser.add_argument('-q', '--quiet', action='store_true',
                            help='quiet mode')
        parser.add_argument('--list', action='store_true',
                            help='print all tests and exit')
        parser.add_argument('--list-calculators', action='store_true',
                            help='print all calculator names and exit')
        parser.add_argument('tests', nargs='*')

    @staticmethod
    def run(args):
        if args.calculators:
            calculators = args.calculators.split(',')
        else:
            calculators = []

        if args.list:
            for testfile in get_tests():
                print(testfile)
            sys.exit(0)

        if args.list_calculators:
            for name in calc_names:
                print(name)
            sys.exit(0)

        for calculator in calculators:
            if calculator not in calc_names:
                sys.stderr.write('No calculator named "{}".\n'
                                 'Possible CALCULATORS are: '
                                 '{}.\n'.format(calculator,
                                                ', '.join(calc_names)))
                sys.exit(1)

        nproblems = test(verbosity=1 + args.verbose - args.quiet,
                         calculators=calculators,
                         files=args.tests)
        sys.exit(nproblems)


if __name__ == '__main__':
    # Run pyflakes on all code in ASE:
    try:
        output = subprocess.check_output(['pyflakes', 'ase', 'doc'],
                                         stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as ex:
        output = ex.output.decode()

    lines = []
    for line in output.splitlines():
        # Ignore these:
        for txt in ['jacapo', 'list comprehension redefines']:
            if txt in line:
                break
        else:
            lines.append(line)
    if lines:
        print('\n'.join(lines))
        sys.exit(1)
