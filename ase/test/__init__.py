from __future__ import print_function
import os
import platform
import sys
import shutil
import subprocess
import tempfile
import unittest
from glob import glob

from ase.calculators.calculator import names as calc_names, get_calculator
from ase.parallel import paropen
from ase.utils import import_module, devnull


class NotAvailable(Exception):
    pass


test_calculator_names = []


def require(calcname):
    if calcname not in test_calculator_names:
        raise NotAvailable
        

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
    def __init__(self, methodname='testfile', filename=None, display=True):
        unittest.TestCase.__init__(self, methodname)
        self.filename = filename
        self.display = display

    def testfile(self):
        try:
            with open(self.filename) as fd:
                exec(compile(fd.read(), self.filename, 'exec'),
                     {'display': self.display})
        except KeyboardInterrupt:
            raise RuntimeError('Keyboard interrupt')
        except ImportError as ex:
            module = ex.args[0].split()[-1].replace("'", '').split('.')[0]
            if module in ['scipy', 'Scientific', 'lxml']:
                sys.__stdout__.write('skipped (no {0} module) '.format(module))
            else:
                raise
        except NotAvailable:
            sys.__stdout__.write('skipped ')

    def id(self):
        return self.filename

    def __str__(self):
        return self.filename.split('test/')[-1]

    def __repr__(self):
        return "ScriptTestCase(filename='%s')" % self.filename


def test(verbosity=1, calculators=[],
         testdir=None, display=True, stream=sys.stdout):
    test_calculator_names.extend(calculators)
    disable_calculators([name for name in calc_names
                         if name not in calculators])
    ts = unittest.TestSuite()
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
    lasttest = None  # is COCu111.py in the current set
    for test in tests:
        if test.endswith('__.py'):
            continue
        if test.endswith('COCu111.py'):
            lasttest = test
            continue
        ts.addTest(ScriptTestCase(filename=os.path.abspath(test),
                                  display=display))
    if lasttest:
        ts.addTest(ScriptTestCase(filename=os.path.abspath(lasttest),
                                  display=display))

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
        print('{0:16}{1}'.format(a, b))
        
    sys.stdout = devnull

    ttr = unittest.TextTestRunner(verbosity=verbosity, stream=stream)

    origcwd = os.getcwd()
    
    if testdir is None:
        testdir = tempfile.mkdtemp(prefix='ase-test-')
    else:
        if os.path.isdir(testdir):
            shutil.rmtree(testdir)  # clean before running tests!
        os.mkdir(testdir)
    os.chdir(testdir)
    print('test-dir       ', testdir, '\n', file=sys.__stdout__)
    try:
        results = ttr.run(ts)
    finally:
        os.chdir(origcwd)
        sys.stdout = sys.__stdout__

    return results


def disable_calculators(names):
    def __init__(self, *args, **kwargs):
        raise NotAvailable

    def __del__(self):
        pass
        
    for name in names:
        if name in ['emt', 'lj', 'eam', 'morse', 'tip3p']:
            continue
        try:
            cls = get_calculator(name)
        except ImportError:
            pass
        else:
            cls.__init__ = __init__
            cls.__del__ = __del__


def cli(command, calculator_name=None):
    if (calculator_name is not None and
        calculator_name not in test_calculator_names):
        return
    error = subprocess.call(' '.join(command.split('\n')), shell=True)
    assert error == 0
    

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

            
if __name__ == '__main__':
    # Run pyflakes3 on all code in ASE:
    try:
        output = subprocess.check_output(['pyflakes3', 'ase', 'doc'])
    except subprocess.CalledProcessError as ex:
        output = ex.output.decode()

    lines = []
    for line in output.splitlines():
        # Ignore these:
        for txt in ['jacapo', 'tasks', 'execute.py',
                    'list comprehension redefines']:
            if txt in line:
                break
        else:
            lines.append(line)
    if lines:
        print('\n'.join(lines))
        sys.exit(1)
