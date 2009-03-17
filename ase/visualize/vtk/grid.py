
import numpy as npy

from vtk import vtkUnstructuredGrid, vtkPoints, vtkIdList, vtkDoubleArray
from cell import vtkUnitCellModule

class vtkAtomicPositions:
    def __init__(self, pos, cell):

        # Make sure position argument is a valid array
        if not isinstance(pos, npy.ndarray):
            pos = npy.array(pos)

        assert pos.dtype == float and pos.shape[1:] == (3,)
        self.positions = pos
        self.npoints = len(pos)

        # Make sure cell argument is correct type
        assert isinstance(cell, vtkUnitCellModule)
        self.cell = cell

        # Transfer atomic positions to VTK points
        self.vtk_pts = vtkPoints()
        self.vtk_pts.SetDataTypeToDouble()
        self.vtk_pts.SetNumberOfPoints(self.npoints)

        for i,pos in enumerate(self.positions):
            self.vtk_pts.InsertPoint(i,pos[0],pos[1],pos[2])

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

    def get_point_data(self):
        return self.vtk_pointdata

    def add_vector_data(self, vec, name=None):

        # Make sure vec argument is a valid array
        if not isinstance(vec, npy.ndarray):
            vec = npy.array(vec)

        assert vec.dtype == float and vec.shape == (self.npoints,3,)

        # Allocate VTK array for vector data
        vtk_vda = vtkDoubleArray()
        vtk_vda.SetNumberOfValues(self.npoints)
        vtk_vda.SetNumberOfComponents(3)
        if name is not None:
            vtk_vda.SetName(name)

        # Transfer vector data to VTK array
        for i,vc in enumerate(vec):
            vtk_vda.InsertTuple3(i,vc[0],vc[1],vc[2])

        # Add VTK array to VTK point data
        self.vtk_pointdata.AddArray(vtk_vda)
        if name is not None:
            self.vtk_pointdata.SetActiveVectors(name)

