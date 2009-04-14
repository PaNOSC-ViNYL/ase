import sys
import unittest
from glob import glob

class NotAvailable(SystemExit):
    def __init__(self, message, code=0):
        SystemExit.__init__(self, (message,code,))
        self.message = message
        self.code = code

if sys.version_info < (2, 4, 0, 'final', 0):
    class CustomTestCase(unittest.TestCase):
        assertTrue = unittest.TestCase.failUnless
        assertFalse = unittest.TestCase.failIf
else:
    from unittest import TestCase as CustomTestCase

class ScriptTestCase(unittest.TestCase):
    def __init__(self, methodname='testfile', filename=None, display=True):
        unittest.TestCase.__init__(self, methodname)
        self.filename = filename
        self.display = display
        
    def testfile(self):
        try:
            execfile(self.filename, {'display': self.display})
        except KeyboardInterrupt:
            raise RuntimeError('Keyboard interrupt')
        except NotAvailable, err:
            # Only non-zero error codes are failures
            if err.code:
                raise err

    def id(self):
        return self.filename

    def __str__(self):
        return '%s (ScriptTestCase)' % self.filename.split('/')[-1]

    def __repr__(self):
        return "ScriptTestCase(filename='%s')" % self.filename


def test(verbosity=1, dir=None, display=True):
    ts = unittest.TestSuite()
    if dir is None:
        dir = __path__[0]
    tests = glob(dir + '/*.py')
    tests.sort()
    for test in tests:
        if test.endswith('__init__.py'):
            continue
        if test.endswith('COCu111.py'):
            lasttest = test
            continue
        ts.addTest(ScriptTestCase(filename=test, display=display))
    ts.addTest(ScriptTestCase(filename=lasttest, display=display))

    from ase.utils import devnull
    sys.stdout = devnull
    
    ttr = unittest.TextTestRunner(verbosity=verbosity, stream=sys.__stdout__)
    results = ttr.run(ts)

    sys.stdout = sys.__stdout__

    return results
