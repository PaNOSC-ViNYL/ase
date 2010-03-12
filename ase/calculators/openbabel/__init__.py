"""
This module contains an interface to the Force Fields available in the
OpenBabel library ( http://openbabel.org/ ).

The OpenBabel Python interface is described in

O'Boyle et al., Chem. Cent. J., 2, 5 (2008), doi:10.1186/1752-153X-2-5
"""

import numpy as np
import openbabel as ob
from ase import units, Atom, Atoms

from ase.calculators.openbabel.tools import atoms_to_obmol

class OBForceField:
    """OpenBabel Force Field calculator

    Currently supporting ghemical and UFF force fields:
    UFF:
        http://openbabel.org/wiki/OBForceFieldUFF
        J. Am. Chem. Soc. 1992, Vol. 114, No. 25, 10024-10035
    ghemical:
        http://openbabel.org/wiki/OBForceFieldGhemical
        http://www.uku.fi/~thassine/projects/ghemical/

    """
    def __init__(self, force_field='UFF', bonds=None):
        """Construct OpenBabel Force Field Calculator object.

        Parameters
        ==========
        force_field: str
            One of 'UFF', 'ghemical'
        bonds: list of lists of 3xint
            Define bonds between atoms such as:
                [[begin atom index, end atom index, bond order],
                 ...
                ]
            If None the calculator will try to construct the bonds
            automatically.
        """

        self.force_field = force_field
        self.atoms = None
        self.bonds = bonds

    def get_potential_energy(self, atoms=None, force_consistent=False):
        """Return total energy"""

        if atoms is None:
            atoms = self.atoms

        mol = atoms_to_obmol(atoms, self.bonds)
        ff = ob.OBForceField.FindForceField(self.force_field)
        ff.Setup(mol)
        E = ff.Energy()

        if ff.GetUnit() == 'kJ/mol':
            E *= units.kJ / units.mol
        elif ff.GetUnit() == 'kcal/mol':
            E *= units.kcal / units.mol
        else:
            raise NotImplementedError

        return E

    def get_forces(self, atoms):
        """Return the forces.

        Forces are currently calculated using finite differences, as openbabel
        does not expose the GetGradient function.
        """

        h = 0.01
        forces = np.zeros_like(atoms.positions)
        for i in range(len(atoms)):
            for j in range(3):
                tmpatoms = atoms.copy()

                tmpatoms.positions[i, j] += h
                forces[i, j] += self.get_potential_energy(tmpatoms)

                tmpatoms.positions[i, j] -= 2 * h
                forces[i, j] -= self.get_potential_energy(tmpatoms)

        # The force is the negative gradient!
        forces /= -(2 * h)

        return forces

    def get_stress(self, atoms):
        """Return the stress."""
        raise NotImplementedError

    def set_atoms(self, atoms):
        self.atoms = atoms.copy()
