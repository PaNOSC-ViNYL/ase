
import numpy as np

from vtk import vtkPolyDataNormals, vtkLinearSubdivisionFilter, \
                vtkSmoothPolyDataFilter, vtkDepthSortPolyData
from ase.visualize.vtk.grid import vtkVolumeGrid

# -------------------------------------------------------------------

class vtkPipeline(list):
    branchable = False

    def __init__(self, vtk_polydata):

        list.__init__(self, [vtk_polydata])

        self.closed = False
        self.pending = False

    def signal_close(self):
        if self.closed:
            raise RuntimeError('Pipeline output port is closed.')
        elif not self.branchable:
            self.pending = True

    def get_output_port(self):
        if self.closed:
            raise RuntimeError('Pipeline output port is closed.')
        elif self.pending:
            self.closed = True
            self.pending = False

        if isinstance(self[-1], vtkPipeline):
            return self[-1].get_output_port()
        else:
            return self[-1].GetOutputPort()

    def connect(self, vtk_polymapper):
        if isinstance(vtk_polymapper, vtkPipeline):
            if not self.branchable and self.closed:
                raise RuntimeError('Pipeline branching not allowed.')

            self.signal_close()
            vtk_polymapper.connect(self.get_output_port())
        else:
            vtk_polymapper.SetInputConnection(self.get_output_port())

    def append(self, vtk_filter):
        self.connect(vtk_filter)
        list.append(self, vtk_filter)

    def extend(self, vtk_filters):
        map(self.append, vtk_filters)

# -------------------------------------------------------------------

class vtkSurfaceSmootherPipeline(vtkPipeline):
    def __init__(self, grid, vtk_polydata, angle=15):

        vtkPipeline.__init__(self, vtk_polydata)

        # Make sure grid argument is correct type
        assert isinstance(grid, vtkVolumeGrid)
        self.grid = grid

        # Split polys with intersection angles greater than angle
        vtk_dnorm = vtkPolyDataNormals()
        vtk_dnorm.SetFeatureAngle(angle)
        vtk_dnorm.SplittingOn()
        vtk_dnorm.ComputeCellNormalsOff()
        vtk_dnorm.ComputePointNormalsOff()
        self.append(vtk_dnorm)

        relax = self.grid.get_relaxation_factor()

        if relax is not None:
            print 'relax=',relax
            #vtk_subd = vtkButterflySubdivisionFilter()
            vtk_subd = vtkLinearSubdivisionFilter()
            self.append(vtk_subd)
    
            # Smooth out some of the sharp points.
            vtk_smooth = vtkSmoothPolyDataFilter()
            vtk_smooth.SetRelaxationFactor(relax)
            self.append(vtk_smooth)

        self.signal_close()

# -------------------------------------------------------------------

class vtkDepthSortPipeline(vtkPipeline):
    def __init__(self, vtk_renderer, vtk_polydata):

        vtkPipeline.__init__(self, vtk_polydata)

        # The depht sort object is set up to generate scalars representing
        # the sort depth. A camera is assigned for the sorting. The camera
        # defines the sort vector (position and focal point).
        vtk_ds = vtkDepthSortPolyData()
        vtk_ds.SetCamera(vtk_renderer.GetActiveCamera())
        vtk_ds.SetDirectionToBackToFront()
        #vtk_ds.SetVector(1, 1, 1)
        #vtk_ds.SortScalarsOn()
        #vtk_ds.Update()
        self.append(vtk_ds)

        vtk_renderer.ResetCamera()

        self.signal_close()

