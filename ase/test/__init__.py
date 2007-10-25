import sys
import unittest
from glob import glob


class ScriptTestCase(unittest.TestCase):
    def __init__(self, methodname='testfile', filename=None):
        unittest.TestCase.__init__(self, methodname)
        self.filename = filename
        
    def testfile(self):
        execfile(self.filename)

    def id(self):
        return self.filename

    def __str__(self):
        return '%s (ScriptTestCase)' % self.filename.split('/')[-1]

    def __repr__(self):
        return "ScriptTestCase(filename='%s')" % self.filename


def test(verbosity=1):
    ts = unittest.TestSuite()
    tests = glob(__path__[0] + '/*.py')
    tests.sort()
    for test in tests:
        if test.endswith('__init__.py'):
            continue
        ts.addTest(ScriptTestCase(filename=test))

    from ase.util import devnull
    sys.stdout = devnull
    
    ttr = unittest.TextTestRunner(verbosity=verbosity)
    ttr.run(ts)

    sys.stdout = sys.__stdout__
