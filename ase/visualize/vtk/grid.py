
import numpy as np

from vtk import vtkPointData, vtkDataArray, vtkUnstructuredGrid, vtkPoints, \
                vtkIdList, vtkStructuredPoints
from ase.visualize.vtk.cell import vtkUnitCellModule
from ase.visualize.vtk.data import vtkDataArrayFromNumPyBuffer, \
                                   vtkDoubleArrayFromNumPyArray, \
                                   vtkDoubleArrayFromNumPyMultiArray

# -------------------------------------------------------------------

class vtkBaseGrid:
    def __init__(self, npoints, cell):
        self.npoints = npoints

        # Make sure cell argument is correct type
        assert isinstance(cell, vtkUnitCellModule)
        self.cell = cell

        self.vtk_pointdata = None

    def set_point_data(self, vtk_pointdata):
        if self.vtk_pointdata is not None:
            raise RuntimeError('VTK point data already present.')

        assert isinstance(vtk_pointdata, vtkPointData)
        self.vtk_pointdata = vtk_pointdata
        self.vtk_pointdata.SetCopyScalars(False)
        self.vtk_pointdata.SetCopyVectors(False)
        self.vtk_pointdata.SetCopyNormals(False)

    def get_point_data(self):
        if self.vtk_pointdata is None:
            raise RuntimeError('VTK point data missing.')

        return self.vtk_pointdata

    def get_number_of_points(self):
        return self.npoints

    def add_scalar_data_array(self, data, name=None, active=True):

        # Are we converting from NumPy buffer to VTK array?
        if isinstance(data, vtkDataArray):
            vtk_sda = data
        elif isinstance(data, vtkDataArrayFromNumPyBuffer):
            vtk_sda = data.get_output()
        else:
            raise ValueError('Data is not a valid scalar data array.')

        del data

        assert vtk_sda.GetNumberOfComponents() == 1
        assert vtk_sda.GetNumberOfTuples() == self.npoints

        if name is not None:
            vtk_sda.SetName(name)

        # Add VTK array to VTK point data
        self.vtk_pointdata.AddArray(vtk_sda)

        if active:
            self.vtk_pointdata.SetActiveScalars(name)

        return vtk_sda

    def add_vector_data_array(self, data, name=None, active=True):

        # Are we converting from NumPy buffer to VTK array?
        if isinstance(data, vtkDataArray):
            vtk_vda = data
        elif isinstance(data, vtkDataArrayFromNumPyBuffer):
            vtk_vda = data.get_output()
        else:
            raise ValueError('Data is not a valid scalar data array.')

        del data

        if name is not None:
            vtk_vda.SetName(name)

        # Add VTK array to VTK point data
        self.vtk_pointdata.AddArray(vtk_vda)

        if active:
            self.vtk_pointdata.SetActiveVectors(name)

        return vtk_vda

# -------------------------------------------------------------------

class vtkAtomicPositions(vtkBaseGrid):
    def __init__(self, pos, cell):

        # Make sure position argument is a valid array
        if not isinstance(pos, np.ndarray):
            pos = np.array(pos)

        assert pos.dtype == float and pos.shape[1:] == (3,)

        vtkBaseGrid.__init__(self, len(pos), cell)

        # Convert positions to VTK array
        npy2da = vtkDoubleArrayFromNumPyArray(pos)
        vtk_pda = npy2da.get_output()
        del npy2da

        # Transfer atomic positions to VTK points
        self.vtk_pts = vtkPoints()
        self.vtk_pts.SetData(vtk_pda)

        # Create a VTK unstructured grid of these points
        self.vtk_ugd = vtkUnstructuredGrid()
        self.vtk_ugd.SetWholeBoundingBox(self.cell.get_bounding_box())
        self.vtk_ugd.SetPoints(self.vtk_pts)

        # Extract the VTK point data set
        self.set_point_data(self.vtk_ugd.GetPointData())

    def get_points(self, subset=None):
        if subset is None:
            return self.vtk_pts

        # Create a list of indices from the subset
        vtk_il = vtkIdList()
        for i in subset:
            vtk_il.InsertNextId(i)

        # Allocate VTK points for subset
        vtk_subpts = vtkPoints()
        vtk_subpts.SetDataType(self.vtk_pts.GetDataType())
        vtk_subpts.SetNumberOfPoints(vtk_il.GetNumberOfIds())

        # Transfer subset of VTK points
        self.vtk_pts.GetPoints(vtk_il, vtk_subpts)

        return vtk_subpts

    def get_unstructured_grid(self, subset=None):
        if subset is None:
            return self.vtk_ugd

        # Get subset of VTK points
        vtk_subpts = self.get_points(subset)

        # Create a VTK unstructured grid of these points
        vtk_subugd = vtkUnstructuredGrid()
        vtk_subugd.SetWholeBoundingBox(self.cell.get_bounding_box())
        vtk_subugd.SetPoints(vtk_subpts)

        return vtk_subugd

    def add_scalar_property(self, data, name=None, active=True):

        # Make sure data argument is a valid array
        if not isinstance(data, np.ndarray):
            data = np.array(data)

        assert data.dtype == float and data.shape == (self.npoints,)

        # Convert scalar properties to VTK array
        npa2da = vtkDoubleArrayFromNumPyArray(data)
        return vtkBaseGrid.add_scalar_data_array(self, npa2da, name, active)

    def add_vector_property(self, data, name=None, active=True):

        # Make sure data argument is a valid array
        if not isinstance(data, np.ndarray):
            data = np.array(data)

        assert data.dtype == float and data.shape == (self.npoints,3,)

        # Convert vector properties to VTK array
        npa2da = vtkDoubleArrayFromNumPyArray(data)
        return vtkBaseGrid.add_vector_data_array(self, npa2da, name, active)

# -------------------------------------------------------------------

class vtkOrthogonalGrid(vtkBaseGrid):
    def __init__(self, elements, cell):

        # Make sure element argument is a valid array
        if not isinstance(elements, np.ndarray):
            elements = np.array(elements)

        assert elements.dtype == int and elements.shape == (3,)
        self.elements = elements

        vtkBaseGrid.__init__(self, np.prod(self.elements), cell)

        # Create a VTK grid of structured points
        self.vtk_spts = vtkStructuredPoints()
        self.vtk_spts.SetWholeBoundingBox(self.cell.get_bounding_box())
        self.vtk_spts.SetDimensions(self.elements)
        self.vtk_spts.SetSpacing(self.get_grid_spacing())

        # Extract the VTK point data set
        self.set_point_data(self.vtk_spts.GetPointData())

    def get_grid_spacing(self):
        return self.cell.get_grid_spacing(self.elements)

    def get_structured_points(self):
        return self.vtk_spts

    def add_scalar_field(self, data, name=None, active=True):

        # Make sure data argument is a valid array
        if not isinstance(data, np.ndarray):
            data = np.array(data)

        assert data.dtype == float and data.shape == tuple(self.elements)

        # Convert scalar field to VTK array
        npa2da = vtkDoubleArrayFromNumPyMultiArray(data[...,np.newaxis])
        return vtkBaseGrid.add_scalar_data_array(self, npa2da, name, active)

    def add_vector_field(self, data, name=None, active=True):

        # Make sure data argument is a valid array
        if not isinstance(data, np.ndarray):
            data = np.array(data)

        assert data.dtype == float and data.shape == tuple(self.elements)+(3,)

        # Convert vector field to VTK array
        npa2da = vtkDoubleArrayFromNumPyMultiArray(data)
        return vtkBaseGrid.add_vector_data_array(self, npa2da, name, active)

