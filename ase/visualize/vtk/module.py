
import numpy as npy

from vtk import vtkProp3D, vtkPolyDataMapper, vtkActor, vtkPointSet, vtkGlyph3D, vtkRenderer
from sources import vtkCustomGlyphSource

class vtkModule:
    def __init__(self, vtk_act, vtk_property=None):

        self.vtk_act = None
        self.set_actor(vtk_act)

        if vtk_property is not None:
            self.set_property(vtk_property)

    def set_actor(self, vtk_act):
        assert isinstance(vtk_act, vtkActor)
        self.vtk_act = vtk_act

    def set_property(self, vtk_property):
        self.vtk_act.SetProperty(vtk_property)

    def get_actor(self):
        return self.vtk_act

# -------------------------------------------------------------------

class vtkPolyDataModule(vtkModule):
    def __init__(self, vtk_polydata, vtk_property=None):

        vtkModule.__init__(self, vtkActor(), vtk_property)

        self.vtk_dmap = vtkPolyDataMapper()
        self.vtk_dmap.SetInputConnection(vtk_polydata.GetOutputPort())
        self.vtk_act.SetMapper(self.vtk_dmap)

class vtkGlyphModule(vtkPolyDataModule):
    def __init__(self, vtk_pointset, glyph_source, clamping=False,
                 scalemode=None, colormode=None):

        assert isinstance(vtk_pointset, vtkPointSet)
        assert isinstance(glyph_source, vtkCustomGlyphSource)

        # Create VTK Glyph3D based on unstructured grid
        self.vtk_g3d = vtkGlyph3D()
        self.vtk_g3d.SetInput(vtk_pointset)
        self.vtk_g3d.SetSource(glyph_source.get_output())
        self.vtk_g3d.SetClamping(clamping)

        if scalemode is 'off':
            self.vtk_g3d.SetScaleModeToDataScalingOff()
        elif scalemode is 'scalar':
            self.vtk_g3d.SetScaleModeToScaleByScalar()
        elif scalemode is 'vector':
            self.vtk_g3d.SetScaleModeToScaleByVector()
        elif scalemode is not None:
            raise ValueError('Inrecognized scale mode \'%s\'.' % scalemode)

        if colormode is 'scale':
            self.vtk_g3d.SetColorModeToColorByScale()
        elif colormode is 'scalar':
            self.vtk_g3d.SetColorModeToColorByScalar()
        elif colormode is 'vector':
            self.vtk_g3d.SetColorModeToColorByVector()
        elif colormode is not None:
            raise ValueError('Inrecognized scale mode \'%s\'.' % scalemode)

        vtkPolyDataModule.__init__(self, self.vtk_g3d, glyph_source.get_property())

# -------------------------------------------------------------------

class vtkLabelModule(vtkModule):
    def __init__(self, vtk_pointset, vtk_property=None):

        vtk_Module.__init__(self, vtkActor(), vtk_property)

        assert isinstance(vtk_pointset, vtkPointSet)

        self.vtk_dmap = vtkLabeledDataMapper()
        #self.vtk_dmap.SetLabelModeToLabelIds() #TODO XXX strings!!!
        self.vtk_dmap.GetLabelTextProperty().SetFontSize(12)
        self.vtk_dmap.GetLabelTextProperty().SetJustificationToRight()
        self.vtk_dmap.SetInputConnection(vtk_pointset.GetOutputPort())
        self.vtk_act.SetMapper(self.vtk_dmap)

        #vtk_g3d.GetPointIdsName...

# -------------------------------------------------------------------

class vtkModuleAnchor:
    def __init__(self):
        self.modules = {}

    def get_module(self, name):
        if name not in self.modules.keys():
            raise RuntimeError('Module \'%s\' does not exists.' % name)

        return self.modules[name]

    def add_module(self, name, module):
        if not isinstance(module, vtkModule):
            raise ValueError('Module must be instance of vtkModule.')

        if name in self.modules.keys():
            raise RuntimeError('Module \'%s\' already exists.' % name)

        self.modules[name] = module

    def get_actor(self, name):
        return self.get_module(name).get_actor()

    def add_actors_to_renderer(self, vtk_renderer, name=None):
        assert isinstance(vtk_renderer, vtkRenderer)
        if name is None:
            for module in self.modules.values():
                vtk_renderer.AddActor(module.get_actor())
        else:
            vtk_renderer.AddActor(self.get_actor(name))

