
import numpy as np

from vtk import vtkContourFilter, vtkDepthSortPolyData
from ase.visualize.vtk.grid import vtkOrthogonalGrid
from ase.visualize.vtk.module import vtkPolyDataModule

class vtkIsoSurfaceModule(vtkOrthogonalGrid, vtkPolyDataModule):
    def __init__(self, data, cell, contours=1, depthsort=True):

        # Make sure data argument is a valid array
        if not isinstance(data, np.ndarray):
            data = np.array(data)

        vtkOrthogonalGrid.__init__(self, data.shape, cell)

        self.vtk_iso = vtkContourFilter() #vtkMarchingContourFilter
        self.vtk_iso.SetInput(self.get_structured_points()) #TODO later, vtkNonOrthogonalGrid common interface here?

        self.vtk_iso.SetValue(0, 0.25) #TODO XXX contour values from list
        self.vtk_iso.SetValue(1, -0.25) #TODO XXX or autocontours if int

        self.depthsort = depthsort

        if self.depthsort:
            # The depht sort object is set up to generate scalars representing
            # the sort depth.  A camera is assigned for the sorting. The camera
            # defines the sort vector (position and focal point).

            self.vtk_ds = vtkDepthSortPolyData()
            self.vtk_ds.SetCamera(ren.GetActiveCamera()) #TODO renderer never passed?!
            self.vtk_ds.SetInputConnection(self.vtk_iso.GetOutputPort())
            self.vtk_ds.SetDirectionToBackToFront()
            #vtk_ds.SetVector(1, 1, 1)
            #vtk_ds.SortScalarsOn()
            #vtk_ds.Update()

            ren.ResetCamera()

            vtkPolyDataModule.__init__(self, self.vtk_ds)
        else:
            vtkPolyDataModule.__init__(self, self.vtk_iso)

