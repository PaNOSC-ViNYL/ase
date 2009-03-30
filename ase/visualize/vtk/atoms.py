"""An experimental package for making plots during a simulation.

A VTKPlotter can plot a list of atoms and most types of volume data.
"""

import numpy as np

from ase import Atoms

from ase.visualize.vtk.sources import vtkAtomSource, vtkForceSource, \
                                      vtkVelocitySource
from ase.visualize.vtk.cell import vtkUnitCellModule, vtkAxesModule
from ase.visualize.vtk.grid import vtkAtomicPositions
from ase.visualize.vtk.module import vtkModuleAnchor, vtkGlyphModule

class vtkAtoms(vtkModuleAnchor, vtkAtomicPositions):
    def __init__(self, atoms, scale=1):

        assert isinstance(atoms, Atoms)
        self.atoms = atoms

        self.scale = scale

        vtkModuleAnchor.__init__(self)
        vtkAtomicPositions.__init__(self, self.atoms.get_positions(),
                                    vtkUnitCellModule(self.atoms))

        self.force = None
        self.velocity = None

        symbols = self.atoms.get_chemical_symbols()
        for symbol in np.unique(symbols):
            # Construct mask for all atoms with this symbol
            mask = np.array(symbols) == symbol
            if mask.all():
                subset = None
            else:
                subset = np.argwhere(mask)

            # Get relevant VTK unstructured grid
            vtk_ugd = self.get_unstructured_grid(subset)

            # Create atomic glyph source for this symbol
            glyph_source = vtkAtomSource(symbol, self.scale)

            # Create glyph module and anchor it
            self.add_module(symbol, vtkGlyphModule(vtk_ugd, glyph_source))

    def has_forces(self):
        return self.force is not None

    def has_velocities(self):
        return self.velocity is not None

    def add_cell(self):
        self.add_module('cell', self.cell)

    def add_axes(self):
        self.add_module('axes', vtkAxesModule(self.cell))

    def add_forces(self):
        if self.has_forces():
            raise RuntimeError('Forces already present.')
        elif self.has_velocities():
            raise NotImplementedError('Can\'t add forces due to velocities.')

        # Add forces to VTK unstructured grid as vector data
        vtk_fda = self.add_vector_property(self.atoms.get_forces(), 'force')

        # Calculate max norm of the forces
        fmax = vtk_fda.GetMaxNorm()

        # Get relevant VTK unstructured grid
        vtk_ugd = self.get_unstructured_grid()

        self.force = vtkGlyphModule(vtk_ugd, vtkForceSource(fmax, self.scale),
                                    scalemode='vector', colormode=None)
        self.add_module('force', self.force)

    def add_velocities(self):
        if self.has_velocities():
            raise RuntimeError('Velocities already present.')
        elif self.has_forces():
            raise NotImplementedError('Can\'t add velocities due to forces.')

        # Add velocities to VTK unstructured grid as vector data
        vtk_vda = self.add_vector_property(self.atoms.get_velocities(), 'velocity')

        # Calculate max norm of the velocities
        vmax = vtk_vda.GetMaxNorm()

        # Get relevant VTK unstructured grid
        vtk_ugd = self.get_unstructured_grid()

        self.velocity = vtkGlyphModule(vtk_ugd, vtkVelocitySource(vmax, self.scale),
                                       scalemode='vector', colormode=None)
        self.add_module('velocity', self.velocity)

