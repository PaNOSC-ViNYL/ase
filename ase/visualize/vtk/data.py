
import numpy as np

from vtk import vtkDataArray, vtkFloatArray, vtkDoubleArray

class vtkNumPyBuffer:
    def __init__(self, data):
        self.strbuf = data.tostring()
        self.nitems = len(data.flat)

    def __len__(self):
        return self.nitems

    def get_pointer(self):
        # Any C/C++ method that requires a void * can be passed a Python
        # string. No check is done to ensure that the string is the correct
        # size, and the string's reference count is not incremented. Extreme
        # caution should be applied when using this feature.
        return self.strbuf

    def notify(self, obj, event):
        if event == 'DeleteEvent':
            del self.strbuf
        else:
            raise RuntimeError('Event not recognized.')

class vtkDataArrayFromNumPyArray:
    """Class for reading vtkDataArray

    This class can be used to generate a vtkDataArray from a NumPy array.
    The NumPy array should be of the form <entries> x <number of components>
    where 'number of components' indicates the number of components in 
    each entry in the vtkDataArray. Note that this form is also expected
    even in the case of only a single component.
    """
    def __init__(self, data, vtk_class, ctype, buffered=True):
        if not isinstance(data, np.ndarray):
            data = np.array(data, dtype=ctype)

        if data.dtype is not ctype:
            data = data.astype(ctype)

        if data.ndim == 1:
            data = data[:, np.newaxis]
        elif data.ndim != 2:
            raise ValueError('Data must be a 2D NumPy array.')

        self.vtk_da = vtk_class()
        assert isinstance(self.vtk_da, vtkDataArray)
        assert self.vtk_da.GetDataTypeSize() == data.itemsize

        self.buffered = buffered

        self.read_numpy_array(data)

    def read_numpy_array(self, data):
        """Read vtkDataArray from NumPy array"""

        if self.buffered:
            self.vtk_da.SetNumberOfComponents(data.shape[-1])

            # Passing the void* buffer to the C interface does not increase
            # its reference count, hence the buffer is deleted by Python when
            # the reference count of the string from tostring reaches zero.
            # Also, the boolean True tells VTK to save (not delete) the buffer
            # when the VTK data array is deleted - we want Python to do this.
            npybuf = vtkNumPyBuffer(data)
            self.vtk_da.SetVoidArray(npybuf.get_pointer(), len(npybuf), True)
            self.vtk_da.AddObserver('DeleteEvent', npybuf.notify)
        else:
            self.vtk_da.SetNumberOfComponents(data.shape[-1])
            self.vtk_da.SetNumberOfTuples(data.shape[0])

            for i, d_c in enumerate(data):
                for c, d in enumerate(d_c):
                    self.vtk_da.SetComponent(i, c, d)

    def get_output(self):
        return self.vtk_da

    def copy(self):
        vtk_da_copy = self.vtk_da.NewInstance()
        vtk_da_copy.SetNumberOfComponents(self.vtk_da.GetNumberOfComponents())
        vtk_da_copy.SetNumberOfTuples(self.vtk_da.GetNumberOfTuples())

        assert vtk_da_copy.GetSize() == self.vtk_da.GetSize()

        vtk_da_copy.DeepCopy(self.vtk_da)

        return vtk_da_copy


class vtkFloatArrayFromNumPyArray(vtkDataArrayFromNumPyArray):
    def __init__(self, data):
        vtkDataArrayFromNumPyArray.__init__(self, data, vtkFloatArray,
                                            np.ctypeslib.ctypes.c_float)

class vtkDoubleArrayFromNumPyArray(vtkDataArrayFromNumPyArray):
    def __init__(self, data):
        vtkDataArrayFromNumPyArray.__init__(self, data, vtkDoubleArray,
                                            np.ctypeslib.ctypes.c_double)

