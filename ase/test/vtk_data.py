#!/usr/bin/env python

from ase.visualize.vtk import requirevtk, probe_vtk_kilobyte

requirevtk()
vtk_kilobyte = probe_vtk_kilobyte(1024)

import unittest
from ase.test import CustomTestCase
import numpy as np

from ase.utils.memory import MemoryStatistics, MemorySingleton, shapeopt

from vtk import vtkDataArray
from ase.visualize.vtk.data import vtkDoubleArrayFromNumPyArray, \
                                   vtkDoubleArrayFromNumPyMultiArray

# Tweak garbage collection (should be on by default)
import gc
gc.enable()

# -------------------------------------------------------------------

class UTConversionDataArrayNumPy(CustomTestCase):
    """
    Abstract class with test cases for VTK/NumPy data conversion.

    Leak tests the six possible permutations of deletion order for the
    objects involved in conversion between VTK and NumPy data formats.
        conv:    The object in charge of the conversion
        data:    NumPy data with a specific memory footprint
        vtk_da:  VTK data array with a similar memory footprint 

    Permutations:
        Case A: 012 i.e. deletion order is conv, data, vtk_da
        Case B: 021 i.e. deletion order is conv, vtk_da, data
        Case C: 102 i.e. deletion order is data, conv, vtk_da
        Case D: 120 i.e. deletion order is data, vtk_da, conv
        Case E: 201 i.e. deletion order is vtk_da, conv, data
        Case F: 210 i.e. deletion order is vtk_da, data, conv
    """
    footprint = 100*1024**2
    verbose = 0
    gc_threshold = (300,5,5) #default is (700,10,10)
    gc_flags = gc.DEBUG_LEAK # | gc.DEBUG_STATS

    def setUp(self):
        self.mem_ini = MemorySingleton(self.verbose-1)
        self.mem_ref = MemoryStatistics(self.verbose-1)
        self.mem_cur = self.mem_ref.copy()

        self.gc_threshold_old = gc.get_threshold()
        self.gc_flags_old = gc.get_debug()
        gc.set_threshold(*self.gc_threshold)
        gc.set_debug(self.gc_flags)

        # Try to obtain a clean slate
        gc.collect()
        self.gc_count = len(gc.garbage)
        del gc.garbage[:]

    def tearDown(self):
        gc.collect()
        self.assertEqual(len(gc.garbage), self.gc_count)
        if len(gc.garbage)>0:
            if self.verbose>1: print gc.get_objects() #DEBUG
            #TODO be pedantic and fail?
        del gc.garbage[:]
        gc.set_threshold(*self.gc_threshold_old)
        gc.set_debug(self.gc_flags_old)

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

    def get_leaktest_scenario(self):
        """Return a tuple of the conversion objects for leak testing.
            conv:    The object in charge of the conversion
            data:    NumPy data with a specific memory footprint
            vtk_da:  VTK data array with a similar memory footprint 
        return (conv, data, vtk_da,)"""

        raise RuntimeError('Virtual member function.')

    # =================================

    def test_deletion_case_a(self):
        # Case A: 012 i.e. deletion order is conv, data, vtk_da
        (conv, data, vtk_da,) = self.get_leaktest_scenario()

        # Conversion cleanup
        del conv
        self.assertAlmostConsumed(2*self.footprint, -5) #100kB -> 140kB sometimes!
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'Conversion cleanup=', self.mem_cur-self.mem_ref #DEBUG

        # NumPy cleanup
        del data
        self.assertAlmostConsumed(self.footprint, -5) #100kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'NumPy cleanup=', self.mem_cur-self.mem_ref #DEBUG

        # VTK cleanup
        del vtk_da
        self.assertAlmostConsumed(0, -4) #10kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'VTK cleanup=', self.mem_cur-self.mem_ref #DEBUG

    def test_deletion_case_b(self):
        # Case B: 021 i.e. deletion order is conv, vtk_da, data
        (conv, data, vtk_da,) = self.get_leaktest_scenario()

        # Conversion cleanup
        del conv
        self.assertAlmostConsumed(2*self.footprint, -5) #100kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'Conversion cleanup=', self.mem_cur-self.mem_ref #DEBUG

        # VTK cleanup
        del vtk_da
        self.assertAlmostConsumed(self.footprint, -5) #100kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'VTK cleanup=', self.mem_cur-self.mem_ref #DEBUG

        # Numpy cleanup
        del data
        self.assertAlmostConsumed(0, -4) #10kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'NumPy cleanup=', self.mem_cur-self.mem_ref #DEBUG

    def test_deletion_case_c(self):
        # Case C: 102 i.e. deletion order is data, conv, vtk_da
        (conv, data, vtk_da,) = self.get_leaktest_scenario()

        # NumPy cleanup
        del data
        self.assertAlmostConsumed(self.footprint, -5) #100kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'NumPy cleanup=', self.mem_cur-self.mem_ref #DEBUG

        # Conversion cleanup
        del conv
        self.assertAlmostConsumed(self.footprint, -5) #100kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'Conversion cleanup=', self.mem_cur-self.mem_ref #DEBUG

        # VTK cleanup
        del vtk_da
        self.assertAlmostConsumed(0, -4) #10kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'VTK cleanup=', self.mem_cur-self.mem_ref #DEBUG

    def test_deletion_case_d(self):
        # Case D: 120 i.e. deletion order is data, vtk_da, conv
        (conv, data, vtk_da,) = self.get_leaktest_scenario()

        # NumPy cleanup
        del data
        self.assertAlmostConsumed(self.footprint, -5) #100kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'NumPy cleanup=', self.mem_cur-self.mem_ref #DEBUG

        # VTK cleanup
        del vtk_da
        self.assertAlmostConsumed(self.footprint, -5) #100kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'VTK cleanup=', self.mem_cur-self.mem_ref #DEBUG

        # Conversion cleanup
        del conv
        self.assertAlmostConsumed(0, -4) #10kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'Conversion cleanup=', self.mem_cur-self.mem_ref #DEBUG

    def test_deletion_case_e(self):
        # Case E: 201 i.e. deletion order is vtk_da, conv, data
        (conv, data, vtk_da,) = self.get_leaktest_scenario()

        # VTK cleanup
        del vtk_da
        self.assertAlmostConsumed(2*self.footprint, -5) #100kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'VTK cleanup=', self.mem_cur-self.mem_ref #DEBUG

        # Conversion cleanup
        del conv
        self.assertAlmostConsumed(self.footprint, -5) #100kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'Conversion cleanup=', self.mem_cur-self.mem_ref #DEBUG

        # NumPy cleanup
        del data
        self.assertAlmostConsumed(0, -4) #10kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'NumPy cleanup=', self.mem_cur-self.mem_ref #DEBUG

    def test_deletion_case_f(self):
        # Case F: 210 i.e. deletion order is vtk_da, data, conv
        (conv, data, vtk_da,) = self.get_leaktest_scenario()

        # VTK cleanup
        del vtk_da
        self.assertAlmostConsumed(2*self.footprint, -5) #100kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'VTK cleanup=', self.mem_cur-self.mem_ref #DEBUG

        # NumPy cleanup
        del data
        self.assertAlmostConsumed(self.footprint, -5) #100kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'NumPy cleanup=', self.mem_cur-self.mem_ref #DEBUG

        # Conversion cleanup
        del conv
        self.assertAlmostConsumed(0, -4) #10kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'Conversion cleanup=', self.mem_cur-self.mem_ref #DEBUG

# -------------------------------------------------------------------

# class UTDataArrayFromNumPyBuffer(...): TODO

# -------------------------------------------------------------------

class UTDataArrayFromNumPyArray_Scalar(UTConversionDataArrayNumPy):
    """
    Test cases for memory leaks during VTK/NumPy data conversion.
    Conversion of 1D NumPy array to VTK data array using buffers."""

    def setUp(self):
        UTConversionDataArrayNumPy.setUp(self)
        self.shape = self.footprint//np.nbytes[np.float]

    def get_leaktest_scenario(self):
        self.assertAlmostEqual(np.prod(self.shape)*np.nbytes[np.float], \
                               self.footprint, -2) #100B

        # Update memory reference
        self.mem_ref.update()

        # NumPy allocation
        data = np.empty(self.shape, np.float)
        self.assertAlmostConsumed(self.footprint, -4) #10kB -> 136kB sometimes!
        self.assertAlmostExceeded(self.footprint, -5) #100kB
        if self.verbose>=1: print 'NumPy allocation=', self.mem_cur-self.mem_ref #DEBUG

        # NumPy to VTK conversion
        np2da = vtkDoubleArrayFromNumPyArray(data)
        self.assertAlmostConsumed(2*self.footprint, -5) #100kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'Conversion=', self.mem_cur-self.mem_ref #DEBUG

        # VTK retrieval
        vtk_da = np2da.get_output()
        self.assertTrue(isinstance(vtk_da, vtkDataArray))
        self.assertAlmostEqual(vtk_da.GetActualMemorySize()*vtk_kilobyte, \
                               self.footprint, -3) #1kB
        if self.verbose>=1: print 'VTK retrieval=', self.mem_cur-self.mem_ref #DEBUG

        return (np2da, data, vtk_da,)

class UTDataArrayFromNumPyArray_Vector(UTConversionDataArrayNumPy):
    """
    Test cases for memory leaks during VTK/NumPy data conversion.
    Conversion of 2D NumPy array to VTK data array using buffers."""

    def setUp(self):
        UTConversionDataArrayNumPy.setUp(self)
        size = self.footprint//np.nbytes[np.float]
        self.shape = (size//3, 3)

    def get_leaktest_scenario(self):
        self.assertAlmostEqual(np.prod(self.shape)*np.nbytes[np.float], \
                               self.footprint, -2) #100B

        # Update memory reference
        self.mem_ref.update()

        # NumPy allocation
        data = np.empty(self.shape, np.float)
        self.assertAlmostConsumed(self.footprint, -4) #10kB
        self.assertAlmostExceeded(self.footprint, -5) #100kB
        if self.verbose>=1: print 'NumPy allocation=', self.mem_cur-self.mem_ref #DEBUG

        # NumPy to VTK conversion
        np2da = vtkDoubleArrayFromNumPyArray(data)
        self.assertAlmostConsumed(2*self.footprint, -5) #100kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'Conversion=', self.mem_cur-self.mem_ref #DEBUG

        # VTK retrieval
        vtk_da = np2da.get_output()
        self.assertTrue(isinstance(vtk_da, vtkDataArray))
        self.assertAlmostEqual(vtk_da.GetActualMemorySize()*vtk_kilobyte, \
                               self.footprint, -3) #1kB
        if self.verbose>=1: print 'VTK retrieval=', self.mem_cur-self.mem_ref #DEBUG

        return (np2da, data, vtk_da,)

# -------------------------------------------------------------------

class UTDataArrayFromNumPyMultiArray_Scalar(UTConversionDataArrayNumPy):
    """
    Test cases for memory leaks during VTK/NumPy data conversion.
    Conversion of NumPy grid scalars to VTK data array using buffers."""

    def setUp(self):
        UTConversionDataArrayNumPy.setUp(self)
        size = self.footprint//np.nbytes[np.float]
        digits, shape = shapeopt(1000, size, ndims=3, ecc=0.3)
        if self.verbose>=1: print 'digits=%8.3f, shape=%s' % (digits,shape) #DEBUG
        self.shape = shape + (1,)
        self.assertAlmostEqual(np.prod(self.shape)*np.nbytes[np.float], \
                               self.footprint, -4) #10kB

    def get_leaktest_scenario(self):

        # Update memory reference
        self.mem_ref.update()

        # NumPy allocation
        data = np.empty(self.shape, np.float)
        self.assertAlmostConsumed(self.footprint, -4) #10kB
        self.assertAlmostExceeded(self.footprint, -5) #100kB
        if self.verbose>=1: print 'NumPy allocation=', self.mem_cur-self.mem_ref #DEBUG

        # NumPy to VTK conversion
        np2da = vtkDoubleArrayFromNumPyMultiArray(data)
        self.assertAlmostConsumed(2*self.footprint, -5) #100kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'Conversion=', self.mem_cur-self.mem_ref #DEBUG

        # VTK retrieval
        vtk_da = np2da.get_output()
        self.assertTrue(isinstance(vtk_da, vtkDataArray))
        self.assertAlmostEqual(vtk_da.GetActualMemorySize()*vtk_kilobyte, \
                               self.footprint, -4) #10kB
        if self.verbose>=1: print 'VTK retrieval=', self.mem_cur-self.mem_ref #DEBUG

        return (np2da, data, vtk_da,)

class UTDataArrayFromNumPyMultiArray_Vector(UTConversionDataArrayNumPy):
    """
    Test cases for memory leaks during VTK/NumPy data conversion.
    Conversion of NumPy grid vectors to VTK data array using buffers."""

    def setUp(self):
        UTConversionDataArrayNumPy.setUp(self)
        size = self.footprint//np.nbytes[np.float]
        digits, shape = shapeopt(1000, size//3, ndims=3, ecc=0.3)
        if self.verbose>=1: print 'digits=%8.3f, shape=%s' % (digits,shape) #DEBUG
        self.shape = shape + (3,)
        self.assertAlmostEqual(np.prod(self.shape)*np.nbytes[np.float], \
                               self.footprint, -4) #10kB

    def get_leaktest_scenario(self):

        # Update memory reference
        self.mem_ref.update()

        # NumPy allocation
        data = np.empty(self.shape, np.float)
        self.assertAlmostConsumed(self.footprint, -4) #10kB
        self.assertAlmostExceeded(self.footprint, -5) #100kB
        if self.verbose>=1: print 'NumPy allocation=', self.mem_cur-self.mem_ref #DEBUG

        # NumPy to VTK conversion
        np2da = vtkDoubleArrayFromNumPyMultiArray(data)
        self.assertAlmostConsumed(2*self.footprint, -5) #100kB
        self.assertAlmostExceeded(2*self.footprint, -6) #1MB
        if self.verbose>=1: print 'Conversion=', self.mem_cur-self.mem_ref #DEBUG

        # VTK retrieval
        vtk_da = np2da.get_output()
        self.assertTrue(isinstance(vtk_da, vtkDataArray))
        self.assertAlmostEqual(vtk_da.GetActualMemorySize()*vtk_kilobyte, \
                               self.footprint, -4) #10kB
        if self.verbose>=1: print 'VTK retrieval=', self.mem_cur-self.mem_ref #DEBUG

        return (np2da, data, vtk_da,)

# -------------------------------------------------------------------

if __name__ in ['__main__', '__builtin__']:
    # We have been imported by test.py, so we should redirect to logfile
    if __name__ == '__builtin__':
        from ase.parallel import paropen
        f = paropen('vtk_data.log', 'w')
    else:
        from sys import stdout as f
    testrunner = unittest.TextTestRunner(verbosity=2, stream=f)

    testcases = [UTDataArrayFromNumPyArray_Scalar, \
                 UTDataArrayFromNumPyArray_Vector, \
                 UTDataArrayFromNumPyMultiArray_Scalar, \
                 UTDataArrayFromNumPyMultiArray_Vector]

    for test in testcases:
        info = '\n' + test.__name__ + '\n' + test.__doc__.strip('\n') + '\n'
        testsuite = unittest.defaultTestLoader.loadTestsFromTestCase(test)
        testrunner.stream.writeln(info)
        testresult = testrunner.run(testsuite)
        # Provide feedback on failed tests if imported by test.py
        if __name__ == '__builtin__' and not testresult.wasSuccessful():
            raise SystemExit('Test failed. Check vtk_data.log for details.')

