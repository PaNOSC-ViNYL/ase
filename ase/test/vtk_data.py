#!/usr/bin/env python

from ase.visualize.vtk import testvtk

testvtk()

import unittest
import numpy as np

from ase.utils.memory import MemoryStatistics, MemorySingleton

from vtk import vtkDataArray
from ase.visualize.vtk.data import vtkDoubleArrayFromNumPyArray, \
                                   vtkDoubleArrayFromNumPyMultiArray

# -------------------------------------------------------------------

class UTNumPyArray(unittest.TestCase):
    """
    Abstract test case class - TODO."""
    footprint = 100*1024**2
    verbose = 1

    def setUp(self):
        self.mem_ini = MemorySingleton(self.verbose-1)
        self.mem_ref = MemoryStatistics(self.verbose-1)
        self.mem_cur = self.mem_ref.copy()

    def assertAlmostConsumed(self, bytes, digits=0, key='VmSize'):
        self.mem_cur.update()
        dm = self.mem_cur-self.mem_ref
        self.assertAlmostEqual(dm[key], bytes, digits)

    def assertAlmostExceeded(self, bytes, digits=0, key='VmPeak'):
        self.mem_cur.update()
        dm = self.mem_cur-self.mem_ini
        #self.assertAlmostEqual(dm[key], bytes, digits) #TODO what really?
        #self.assertAlmostEqual(max(0, dm[key]-bytes), 0, digits) #TODO ???
        #dm = 200 MB, bytes = 100MB     ok
        #dm = 101 MB, bytes = 100MB     ok
        #dm = 0 MB, bytes = 100MB       bad

    def test_consistency_simple(self):
        size = self.footprint//np.nbytes[np.float]
        self.assertEqual(size*np.nbytes[np.float], self.footprint, -1) #10B

        # Reference
        self.mem_ref.update()

        # NumPy allocation
        data = np.empty(size, np.float)
        self.assertAlmostConsumed(self.footprint, -4) #10kB
        self.assertAlmostExceeded(self.footprint, -5) #100kB
        if self.verbose>=1: print 'NumPy allocation=', self.mem_cur-self.mem_ref #DEBUG

        # NumPy to VTK conversion
        np2da = vtkDoubleArrayFromNumPyArray(data)
        vtk_da = np2da.get_output()
        del np2da
        self.assertTrue(isinstance(vtk_da, vtkDataArray))
        self.assertAlmostEqual(vtk_da.GetActualMemorySize()*1e3, \
                               self.footprint, -3) #1kB
        self.assertAlmostConsumed(2*self.footprint, -5) #100kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'VTK allocation=', self.mem_cur-self.mem_ref #DEBUG

        # Numpy cleanup
        del data
        self.assertAlmostConsumed(self.footprint, -5) #100kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'NumPy cleanup=', self.mem_cur-self.mem_ref #DEBUG

        # VTK cleanup
        del vtk_da
        self.assertAlmostConsumed(0, -4) #10kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'VTK cleanup=', self.mem_cur-self.mem_ref #DEBUG

# -------------------------------------------------------------------

if __name__ == '__main__':
    testrunner = unittest.TextTestRunner(verbosity=2)

    testcases = [UTNumPyArray]

    for test in testcases:
        info = '\n' + test.__name__ + '\n' + test.__doc__.strip('\n') + '\n'
        testsuite = unittest.defaultTestLoader.loadTestsFromTestCase(test)
        testrunner.stream.writeln(info)
        testrunner.run(testsuite)

