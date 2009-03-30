
import numpy as np

from ase import Atoms

from vtk import vtkOutlineSource, vtkAxesActor, vtkProperty2D, vtkTextProperty
from ase.visualize.vtk.module import vtkModule, vtkPolyDataModule

class vtkUnitCellModule(vtkPolyDataModule):
    def __init__(self, atoms):

        assert isinstance(atoms, Atoms)
        cell = atoms.get_cell()

        """
        if not isinstance(cell, np.ndarray):
            cell = np.array(cell)

        if cell.shape == (3,):
            cell = np.diag(cell)

        assert cell.dtype == float and cell.shape == (3, 3)
        """

        #TODO bounding box with general unit cell?!
        diagcell = np.diag(cell.diagonal())
        assert (cell == diagcell).all(), 'Unit cell must be orthogonal'

        self.bbox = np.array(zip(np.zeros(3),cell.diagonal())).ravel()

        self.pbc = atoms.get_pbc()

        self.vtk_outline = vtkOutlineSource()
        self.vtk_outline.SetBounds(self.bbox)

        vtkPolyDataModule.__init__(self, self.vtk_outline)

    def get_bounding_box(self):
        return self.bbox

    def get_size(self):
        return self.bbox[1::2]-self.bbox[0::2]

    def get_characteristic_length(self):
        return np.prod(self.get_size())**(1.0/3.0)

    def get_grid_spacing(self, shape):
        return self.get_size()/(np.array(shape)-1.0) #TODO pbc

    def get_relaxation_factor(self, shape):
        # The relaxation factor is a floating point value between zero and one.
        # It expresses the need for smoothening (relaxation) e.g. of isosurfaces
        # due to coarse grid spacings. Larger grid spacing -> larger relaxation.
        x = self.get_grid_spacing(shape).mean()/self.get_characteristic_length()

        # The relaxation function f(x) satisfies the following requirements
        # f(x) -> 0 for x -> 0+   and   f(x) -> b for x -> inf
        # f'(x) -> a for x -> 0+  and   f'(x) -> 0 for x -> inf

        # Furthermore, it is a rescaling of arctan, hence we know
        # f(x) = 2 b arctan(a pi x / 2 b) / pi

        # Our reference point is x = r for which medium relaxion is needed
        # f(r) = b/2   <=>   r = 2 b / a pi   <=>   a = 2 b / r pi
        r = 0.025 # corresponding to 0.2 Ang grid spacing in 8 Ang cell
        b = 0.5
        f = 2*b*np.arctan(x/r)/np.pi

        if f > 0.1:
           return f.round(1)
        else:
           return None

# -------------------------------------------------------------------

class vtkAxesModule(vtkModule):
    def __init__(self, cell):

        assert isinstance(cell, vtkUnitCellModule)
        self.cell = cell

        l0 = self.cell.get_characteristic_length()

        # Create VTK axes actor (not really a VTK actor though)
        self.vtk_ax = vtkAxesActor()
        self.vtk_ax.SetTipTypeToCone()
        self.vtk_ax.SetConeRadius(5e-2*l0)
        self.vtk_ax.SetShaftTypeToCylinder()
        self.vtk_ax.SetCylinderRadius(5e-3*l0)

        # Create VTK two-dimensional property
        p2d = vtkProperty2D()
        p2d.SetDisplayLocationToBackground()

        vtkModule.__init__(self, self.vtk_ax, p2d)

        # Create VTK text property and apply to axes
        vtk_textproperty = vtkTextProperty()
        vtk_textproperty.SetFontSize(14)
        vtk_textproperty.SetBold(True)
        vtk_textproperty.SetItalic(True)
        vtk_textproperty.SetShadow(True)
        vtk_textproperty.SetJustificationToRight()
        vtk_textproperty.SetVerticalJustificationToCentered()

        self.set_text_property(vtk_textproperty)

    def set_actor(self, vtk_act):
        assert isinstance(vtk_act, vtkAxesActor) #fix for non-vtkActor actor
        self.vtk_act = vtk_act

    def set_property(self, vtk_property):
        assert isinstance(vtk_property, vtkProperty2D)
        for vtk_cap in [self.vtk_ax.GetXAxisCaptionActor2D(),
                        self.vtk_ax.GetYAxisCaptionActor2D(),
                        self.vtk_ax.GetZAxisCaptionActor2D()]:
            #vtk_cap.ThreeDimensionalLeaderOn()
            #vtk_cap.LeaderOn()
            vtk_cap.SetProperty(vtk_property)
            vtk_txt = vtk_cap.GetTextActor()
            vtk_txt.SetProperty(vtk_property)

    def set_text_property(self, vtk_textproperty, scaled=False):
        assert isinstance(vtk_textproperty, vtkTextProperty)
        for vtk_cap in [self.vtk_ax.GetXAxisCaptionActor2D(),
                        self.vtk_ax.GetYAxisCaptionActor2D(),
                        self.vtk_ax.GetZAxisCaptionActor2D()]:
            vtk_txt = vtk_cap.GetTextActor()
            vtk_txt.SetScaledText(scaled)
            vtk_txt.SetTextProperty(vtk_textproperty)

