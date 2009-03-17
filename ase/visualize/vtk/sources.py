
import numpy as npy

from vtk import vtkProperty, vtkSphereSource, vtkArrowSource, vtkConeSource

from ase.data import atomic_numbers
from ase.data import covalent_radii as atomic_radii
from ase.data.colors import jmol_colors as atomic_colors
#from ase.data.colors import cpk_colors as atomic_colors


class vtkCustomGlyphSource:
    def __init__(self, scale, glyph_source):
        self.scale = scale
        self.vtk_glyph_source = glyph_source
        self.vtk_property = vtkProperty()

    def get_property(self):
        return self.vtk_property

    def get_output(self):
        return self.vtk_glyph_source.GetOutput()

class vtkAtomSource(vtkCustomGlyphSource):
    def __init__(self, name, scale=1):
        vtkCustomGlyphSource.__init__(self, scale, vtkSphereSource())

        self.number = atomic_numbers[name]
        self.radius = atomic_radii[self.number]
        self.color = atomic_colors[self.number]

        self.vtk_property.SetColor(self.color[0],self.color[1],self.color[2])
        self.vtk_property.SetInterpolationToPhong()
        self.vtk_property.SetDiffuse(0.7)
        self.vtk_property.SetSpecular(0.4)
        self.vtk_property.SetSpecularPower(20)

        self.vtk_glyph_source.SetPhiResolution(16)
        self.vtk_glyph_source.SetThetaResolution(16)
        self.vtk_glyph_source.SetRadius(self.scale*self.radius)

class vtkForceSource(vtkCustomGlyphSource):
    def __init__(self, scale=1):
        vtkCustomGlyphSource.__init__(self, scale, vtkArrowSource())

        self.vtk_property.SetInterpolationToPhong()
        self.vtk_property.SetDiffuse(0.7)
        self.vtk_property.SetSpecular(0.4)
        self.vtk_property.SetSpecularPower(20)

        self.vtk_glyph_source.SetShaftResolution(12)
        self.vtk_glyph_source.SetTipResolution(20)

class vtkVelocitySource(vtkCustomGlyphSource):
    def __init__(self, scale=1):
        vtkCustomGlyphSource.__init__(self, scale, vtkConeSource())

        self.vtk_glyph_source.SetAngle(10)
        self.vtk_glyph_source.SetResolution(16)

class vtkBondSource(vtkCustomGlyphSource):
    def __init__(self, width, scale=1):
        vtkCustomGlyphSource.__init__(self, scale, vtkCylinderSource())

        self.width = width

        self.vtk_glyph_source.SetRadius(self.scale*self.width)
        self.vtk_glyph_source.SetResolution(16)

