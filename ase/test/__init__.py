from __future__ import print_function
import os
import sys
import subprocess
from multiprocessing import Process, cpu_count, Queue
try:
    import queue
except ImportError:
    import Queue as queue  # Python2
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


def runtest_almost_no_magic(test):
    try:
        with open(test) as fd:
            exec(compile(fd.read(), test, 'exec'), {})
    except ImportError as ex:
        module = ex.args[0].split()[-1].replace("'", '').split('.')[0]
        if module in ['scipy', 'matplotlib', 'Scientific', 'lxml',
                      'flask', 'gpaw', 'GPAW', 'netCDF4']:
            raise unittest.SkipTest('no {} module'.format(module))
        else:
            raise


def run_single_test(filename):
    """Execute single test and return results as dictionary."""
    sys.stdout = devnull
    basedir = os.path.split(__file__)[0]
    testrelpath = os.path.relpath(filename, basedir)
    result = {'name': testrelpath, 'pid': os.getpid()}
    #test = ScriptTestCase(filename=filename)

    # Some tests may write to files with the same name as other tests.
    # Hence, create new subdir for each test:
    cwd = os.getcwd()
    testsubdir = testrelpath.replace('/', '_').replace('.', '_')
    os.mkdir(testsubdir)
    os.chdir(testsubdir)
    t1 = time.time()

    try:
        runtest_almost_no_magic(filename)
    except KeyboardInterrupt:
        raise
    except unittest.SkipTest as ex:
        result['status'] = 'SKIPPED'
        result['whyskipped'] = str(ex)
        result['ex'] = ex
    except AssertionError as ex:
        result['status'] = 'FAIL'
        result['ex'] = ex
        result['traceback'] = traceback.format_exc()
    except BaseException as ex:  # Return any error/signal to parent process.
        # Important: In Py2, subprocesses get the Ctrl+C signal whereas in
        # Py3, the parent process does - which is much better.
        # In Py2 it is difficult to kill the processes since the parent
        # process does not know about Ctrl+C.  This is why we
        # return the exception to the parent process and let the
        # parent process smoothly close the pool.
        result['status'] = 'ERROR'
        result['ex'] = ex
        result['traceback'] = traceback.format_exc()
    else:
        result['status'] = 'OK'
    finally:
        sys.stdout = sys.__stdout__
        t2 = time.time()
        os.chdir(cwd)

    result['time'] = t2 - t1
    return result


def runtests_subprocess(task_queue, result_queue):
    """Main test loop to be called within subprocess."""
    test = 'N/A'

    try:
        while True:
            result = None
            try:
                test = task_queue.get_nowait()
            except queue.Empty:
                return  # Done!

            if 'run/gui.py' in test:
                result_queue.put('PLEASE RUN THIS TEST ON MASTER')
                continue

            result = run_single_test(test)
            result_queue.put(result)
    except KeyboardInterrupt:
        result = None
        return
    except BaseException as err:
        # Failure outside actual test (so test suite error):
        result = {'pid': os.getpid(), 'name': test, 'ex': err,
                  'traceback': traceback.format_exc(),
                  'time': 0.0, 'status': 'TESTSUITEBROKEN'}
    finally:
        pass  # finalize things somehow
        #if result is None:
        #    result = 'bogus result'
        #result_queue.put(result)


class Test:
    testbasedir = os.path.split(__file__)[0]

    def __init__(self, name):
        self.fullpath = name
        self.name = os.path.relpath(self.fullname, basedir)


def print_test_result(result, i, n, tests, results):
    msg = result['status']
    if msg == 'SKIPPED':
        msg = 'SKIPPED: {}'.format(result['whyskipped'])
    finished_tests = set([r['name'] for r in results])
    basedir = os.path.split(__file__)[0]
    reltests = [os.path.relpath(t, basedir) for t in tests]

    stillnotfinished = [test for test in reltests if test not in finished_tests]
    print('{pid:5d}: {name:28} {time:6.2f}s {msg} ({i}/{n}) [remain: {N} ({X})]'
          .format(msg=msg, i=i, n=n, N=len(stillnotfinished),
                  X=','.join(stillnotfinished[:4]), **result))
    if 'traceback' in result:
        print('=' * 78)
        print('Error in {} from pid {}:'.format(result['name'], result['pid']))
        print(result['traceback'].rstrip())
        print('=' * 78)


def runtests_parallel(nprocs, tests):
    # Put all tests in a synchronized queue:
    task_queue = Queue()

    for test in tests:
        # gui/run won't work in other processes.  Nobody knows why
        task_queue.put(test)


    # Test results to be received into other queue:
    result_queue = Queue()
    results = []

    procs = []
    try:
        # Start tasks:
        for i in range(nprocs):
            p = Process(target=runtests_subprocess,
                        name='ASE-test-worker-{}'.format(i),
                        args=[task_queue, result_queue])
            procs.append(p)
            p.start()

        # Collect results:
        for i in range(len(tests)):
            result = result_queue.get()  # blocking call
            if result == 'PLEASE RUN THIS TEST ON MASTER':
                result = run_single_test(test)
            print_test_result(result, i, len(tests), tests, results)
            results.append(result)
            #if result['status'] == 'TESTSUITEBROKEN':
            #    raise result['ex']
        # TODO: Print summary here
    except BaseException as err:
        for proc in procs:
            proc.terminate()
            proc.join()
        del procs[:]
        traceback.print_exc()
    else:
        #for test in tests_master_only:
        #    result = run_single_test(test)
        #    results.append(result)
        return results
    finally:
        for proc in procs:
            proc.join()



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
    nprocs = cpu_count()

    print_info()

    origcwd = os.getcwd()
    testdir = tempfile.mkdtemp(prefix='ase-test-')
    os.chdir(testdir)

    # Note: :25 corresponds to ase.cli indentation
    print('{:25}{}'.format('test directory', testdir))
    print('{:25}{}'.format('number of processes', nprocs))
    print('{:25}{}'.format('time', time.strftime('%c')))
    print()

    t1 = time.time()
    results = []
    try:
        for result in runtests_parallel(nprocs, tests):
            results.append(result)
    except KeyboardInterrupt:
        print('\nInterrupted by keyboard')
    finally:
        t2 = time.time()
        print('Time elapsed: {:.1f} s'.format(t2 - t1))
        trouble = [r for r in results if r['status'] in ['FAIL', 'ERROR']]
        return trouble


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

        trouble = test(verbosity=1 + args.verbose - args.quiet,
                       calculators=calculators,
                       files=args.tests)
        sys.exit(len(trouble))


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
