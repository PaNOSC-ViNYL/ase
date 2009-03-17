"""An experimental package for making plots during a simulation.

A VTKPlotter can plot a list of atoms and most types of volume data.
"""

import numpy as npy

from ase import Atoms

from sources import vtkAtomSource, vtkForceSource, vtkVelocitySource
from cell import vtkUnitCellModule, vtkAxesModule
from grid import vtkAtomicPositions
from module import vtkModuleAnchor, vtkGlyphModule

class vtkAtoms(vtkModuleAnchor, vtkAtomicPositions):
    def __init__(self, atoms, scale=0.25):

        assert isinstance(atoms, Atoms)
        self.atoms = atoms

        vtkModuleAnchor.__init__(self)
        vtkAtomicPositions.__init__(self, self.atoms.get_positions(),
                                    vtkUnitCellModule(self.atoms))

        self.forces = None
        self.velocities = None

        symbols = self.atoms.get_chemical_symbols()
        for symbol in npy.unique(symbols):
            # Construct mask for all atoms with this symbol
            mask = npy.array(symbols) == symbol
            if mask.all():
                subset = None
            else:
                subset = npy.argwhere(mask)

            # Get relevant VTK unstructured grid
            vtk_ugd = self.get_unstructured_grid(subset)

            # Create atomic glyph source for this symbol
            glyph_source = vtkAtomSource(symbol, scale)

            # Create glyph module and anchor it
            self.add_module(symbol, vtkGlyphModule(vtk_ugd, glyph_source))

    def has_forces(self):
        return self.forces is not None

    def has_velocities(self):
        return self.velocities is not None

    """
    def get_glyph_source(self, symbol):
        return self.glyph_sources[symbol]

    def get_point_collection(self, symbol):
        return self.point_collections[symbol]
    """

    def add_cell(self):
        self.add_module('cell', self.cell)

    def add_axes(self):
        self.add_module('axes', vtkAxesModule(self.cell))

    def add_forces(self):
        if self.has_forces():
            raise RuntimeError('Forces already present.')

        # Add forces to VTK unstructured grid as vector data
        self.add_vector_data(self.atoms.get_forces(), 'forces')

        # Get relevant VTK unstructured grid
        vtk_ugd = self.get_unstructured_grid()

        self.forces = vtkGlyphModule(vtk_ugd, vtkForceSource(), clamping=True,
                                     scalemode='vector', colormode=None)
        self.add_module('forces', self.forces)

    def add_velocities(self):
        if self.has_velocities():
            raise RuntimeError('Velocities already present.')

        # Add forces to VTK unstructured grid as vector data
        self.add_vector_data(self.atoms.get_velocities(), 'velocities')

        # Get relevant VTK unstructured grid
        vtk_ugd = self.get_unstructured_grid()

        self.velocities = vtkGlyphModule(vtk_ugd, vtkVelocitySource(), clamping=True,
                                         scalemode='vector', colormode=None)
        self.add_module('velocities', self.velocities) #TODO XXX active vector clash!

