
import numpy as np

from vtk import vtkPointData, vtkUnstructuredGrid, vtkPoints, vtkIdList, \
                vtkStructuredPoints
from ase.visualize.vtk.cell import vtkUnitCellModule
from ase.visualize.vtk.data import vtkDoubleArrayFromNumPyArray

# -------------------------------------------------------------------

class vtk3DGrid:
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

    def get_point_data(self):
        if self.vtk_pointdata is None:
            raise RuntimeError('VTK point data missing.')

        return self.vtk_pointdata

    def get_number_of_points(self):
        return self.npoints

    def add_scalar_data(self, sca, name=None):

        # Make sure sca argument is a valid array
        if not isinstance(sca, np.ndarray):
            sca = np.array(sca)

        assert sca.dtype == float and sca.shape == (self.npoints,)

        # Convert positions to VTK array
        npy2da = vtkDoubleArrayFromNumPyArray(sca)
        vtk_sda = npy2da.get_output()
        del npy2da

        if name is not None:
            vtk_sda.SetName(name)

        # Add VTK array to VTK point data
        self.vtk_pointdata.AddArray(vtk_sda)
        if name is not None:
            self.vtk_pointdata.SetActiveScalars(name)

        return vtk_sda

    def add_vector_data(self, vec, name=None):

        # Make sure vec argument is a valid array
        if not isinstance(vec, np.ndarray):
            vec = np.array(vec)

        assert vec.dtype == float and vec.shape == (self.npoints,3,)

        # Convert positions to VTK array
        npy2da = vtkDoubleArrayFromNumPyArray(vec)
        vtk_vda = npy2da.get_output()
        del npy2da

        if name is not None:
            vtk_vda.SetName(name)

        # Add VTK array to VTK point data
        self.vtk_pointdata.AddArray(vtk_vda)
        if name is not None:
            self.vtk_pointdata.SetActiveVectors(name)

        return vtk_vda

# -------------------------------------------------------------------

class vtkAtomicPositions(vtk3DGrid):
    def __init__(self, pos, cell):

        # Make sure position argument is a valid array
        if not isinstance(pos, np.ndarray):
            pos = np.array(pos)

        assert pos.dtype == float and pos.shape[1:] == (3,)

        vtk3DGrid.__init__(self, len(pos), cell)

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

        # Create a VTK point data set
        self.vtk_pointdata = self.vtk_ugd.GetPointData()

    def get_points(self, subset=None):
        if subset is None:
            return self.vtk_pts

        # Create a list of indices from the subset
        vtk_il = vtkIdList()
        for i in subset:
            vtk_il.InsertNextId(i)

        # Allocate VTK points for subset
        vtk_subpts = vtkPoints()
        vtk_subpts.SetDataTypeToDouble()
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

# -------------------------------------------------------------------

class vtkOrthogonalGrid(vtk3DGrid):
    def __init__(self, elements, cell):

        # Make sure element argument is a valid array
        if not isinstance(elements, np.ndarray):
            elements = np.array(elements)

        assert elements.dtype == int and elements.shape == (3,)
        self.elements = elements

        vtk3DGrid.__init__(self, np.prod(self.elements), cell)

        # Create a VTK grid of structured points
        self.vtk_spts = vtkStructuredPoints()
        self.vtk_spts.SetWholeBoundingBox(self.cell.get_bounding_box())
        self.vtk_spts.SetDimensions(self.elements)
        self.vtk_spts.SetSpacing(self.get_grid_spacing())

        # Assign the VTK array as point data of the grid
        self.vtk_pointdata = self.vtk_spts.GetPointData()

    def get_grid_spacing(self):
        return self.cell.get_grid_spacing(self.elements)

    def get_structured_points(self):
        return self.vtk_spts

    def add_scalar_data(self, sca, name=None):
        # Make sure sca argument is a valid array
        if not isinstance(sca, np.ndarray):
            sca = np.array(sca)

        assert sca.dtype == float and sca.shape == tuple(self.elements)

        vtk3DGrid.add_scalar_data(self, sca.swapaxes(0,2).ravel(), name)

