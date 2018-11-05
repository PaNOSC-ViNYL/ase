"""
Provides FixSymmetry class to preserve spagegroup symmetry during optimisation
"""
from __future__ import print_function
import sys
import numpy as np

from ase.constraints import (FixConstraint,
                             voigt_6_to_full_3x3_stress,
                             full_3x3_to_voigt_6_stress)

# check if we have access to get_spacegroup from spglib
# https://atztogo.github.io/spglib/
has_spglib = False
try:
    import spglib                   # For version 1.9 or later
    has_spglib = True
except ImportError:
    try:
        from pyspglib import spglib  # For versions 1.8.x or before
        has_spglib = True
    except ImportError:
        pass

__all__ = ['refine', 'check', 'FixSymmetry']

def refine(at, symprec=0.01):
    # test orig config with desired tol
    dataset = spglib.get_symmetry_dataset(at, symprec=symprec)
    if dataset is None:
        raise ValueError("refine failed to get initial symmetry dataset " +
                         spglib.get_error_message())
    print("symmetry.refine_symmetry: loose ({}) initial symmetry group number"
          "{}, international (Hermann-Mauguin) {} Hall {}".format(symprec,
          dataset["number"],dataset["international"],dataset["hall"])

    # make sure to use a consistent set of primitive cell,
    # mapping to primitive, and (in the next section)
    # rotation matrix derived from standard cell and transformation matrix
    (prim_cell, prim_scaled_pos, prim_types) = spglib.find_primitive(at,
                                                               symprec=symprec)
    mapping = dataset["mapping_to_primitive"]

    # symmetrize cell
    # create symmetrized rotated cell from standard lattice
    # and transformation matrix
    transformed_std_cell = np.dot(dataset['std_lattice'].T,
                                  dataset['transformation_matrix']).T
    # SVD decompose transformation from symmetrized cell to original
    # cell to extract rotation
    # see https://igl.ethz.ch/projects/ARAP/svd_rot.pdf Eq. 17 and 21
    # S = X Y^T, X = symm_cell.T, Y = orig_cell.T
    (u, s, v_T) = np.linalg.svd(np.dot(transformed_std_cell.T, at.get_cell()))
    rotation = np.dot(v_T.T, u.T)
    # align symmetrized cell by rotating transformed standard cell
    aligned_transformed_std_cell = np.dot(transformed_std_cell, rotation.T)
    # set new cell to aligned symmetrized cell (approx. aligned with
    # original directions, not axes)
    at.set_cell(aligned_transformed_std_cell, True)

    #rotate primitive cell (the same rotation as for the standard cell)
    # to align with orig cell
    aligned_prim_cell = np.dot(prim_cell, rotation.T)
    aligned_prim_pos = np.dot(prim_scaled_pos, aligned_prim_cell)

    # calculate inverse prim cell matrix for transformations
    aligned_prim_cell_inv = np.linalg.inv(aligned_prim_cell)

    # vector to align atom 0 with corresponding atom in primitive cell
    pos = at.get_positions()
    alignment_vec = pos[0] - aligned_prim_pos[mapping[0]]

    symm_pos = pos.copy()
    for i_at in range(len(at)):
        # pos =~ L . prim_cell + prim_pos[mapping] + alignment_vec
        # where L is a vector of integers
        # L = round( (pos - alignment_vec - prim_pos[mapping]) . prim_cell^-1 )
        L_raw = np.dot(pos[i_at] - alignment_vec -
                       aligned_prim_pos[mapping[i_at]], aligned_prim_cell_inv)
        L = np.round(L_raw)
        if np.linalg.norm(np.dot(L-L_raw,aligned_prim_cell)) > symprec*2.0:
            raise ValueError("non-integer distance (in lattice coords) between "
                             "original pos and corresponding primitive"
                             " cell {}   {} {} {}".
                             format(i_at,L_raw[0],L_raw[1],L_raw[2]))
        symm_pos[i_at] = (np.dot(L, aligned_prim_cell) +
                                aligned_prim_pos[mapping[i_at]] +
                                alignment_vec)

    at.set_positions(symm_pos)

    # test final config with tight tol
    dataset = spglib.get_symmetry_dataset(at, symprec=1.0e-4)
    if dataset is None:
        raise ValueError("refine failed to get final symmetry dataset " +
                         spglib.get_error_message())
    print("symmetry.refine_symmetry: precise ({}) symmetrized symmetry group "
          "number {}, international (Hermann-Mauguin) {} Hall {}".format(1.0e-4,
          dataset["number"],dataset["international"],dataset["hall"])

def check(at, symprec=1.0e-6):
    """
    Check symmetry of `at` with precision `symprec` using `spglib`

    Prints a summary and returns result of `spglib.get_symmetry_dataset()`
    """
    dataset = spglib.get_symmetry_dataset(at, symprec=symprec)
    print("ase.spacegroup.symmetrize.check: prec", symprec,
          "got symmetry group number", dataset["number"],
          ", international (Hermann-Mauguin)", dataset["international"],
          ", Hall ",dataset["hall"])
    return dataset

def prep(at, symprec=1.0e-6):
    """
    Prepare `at` for symmetry-preserving minimisation at precision `symprec`

    Returns a tuple `(rotations, translations, symm_map)`
    """
    dataset = spglib.get_symmetry_dataset(at, symprec=symprec)
    print("symmetry.prep: symmetry group number",dataset["number"],
          ", international (Hermann-Mauguin)", dataset["international"],
          ", Hall", dataset["hall"])
    rotations = dataset['rotations'].copy()
    translations = dataset['translations'].copy()
    symm_map=[]
    scaled_pos = at.get_scaled_positions()
    for (r, t) in zip(rotations, translations):
        this_op_map = [-1] * len(at)
        for i_at in range(len(at)):
            new_p = np.dot(r, scaled_pos[i_at,:]) + t
            dp = scaled_pos - new_p
            dp -= np.round(dp)
            i_at_map = np.argmin(np.linalg.norm(dp,  axis=1))
            this_op_map[i_at] = i_at_map
        symm_map.append(this_op_map)
    return (rotations, translations, symm_map)

def symmetrize_forces(lattice, inv_lattice, forces, rot, trans, symm_map):
    """
    Return symmetrized forces

    lattice vectors expected as row vectors (same as ASE get_cell() convention),
    inv_lattice is its matrix inverse (get_reciprocal_cell().T)
    """
    scaled_symmetrized_forces_T = np.zeros(forces.T.shape)

    scaled_forces_T = np.dot(inv_lattice.T,forces.T)
    for (r, t, this_op_map) in zip(rot, trans, symm_map):
        transformed_forces_T = np.dot(r, scaled_forces_T)
        scaled_symmetrized_forces_T[:,this_op_map[:]] +=
            transformed_forces_T[:,:]
    scaled_symmetrized_forces_T /= len(rot)

    symmetrized_forces = np.dot(lattice.T, scaled_symmetrized_forces_T).T

    return symmetrized_forces

def symmetrize_stress(lattice, lattice_inv, stress_3_3, rot):
    """
    Return symmetrized stress

    lattice vectors expected as row vectors (same as ASE get_cell() convention),
    inv_lattice is its matrix inverse (get_reciprocal_cell().T)
    """
    scaled_stress = np.dot(np.dot(lattice, stress_3_3), lattice.T)

    symmetrized_scaled_stress = np.zeros((3,3))
    for r in rot:
        symmetrized_scaled_stress += np.dot(np.dot(r.T, scaled_stress),r)
    symmetrized_scaled_stress /= len(rot)

    return np.dot(np.dot(lattice_inv, symmetrized_scaled_stress),
                  lattice_inv.T)

class FixSymmetry(FixConstraint):
    """
    Constraint to preserve spacegroup symmetry during optimisation.

    Requires spglib package to be available.
    """
    def __init__(self, atoms):
        if not has_spglib:
            import spglib # generate ImportError to skip tests
        self.rotations, self.translations, self.symm_map = prep(atoms)

    def adjust_cell(self, atoms, cell):
        pass

    def adjust_positions(self, atoms, new):
        pass

    def adjust_forces(self, atoms, forces):
        forces[:] = symmetrize_forces(atoms.get_cell(),
                                      atoms.get_reciprocal_cell().T,
                                      forces,
                                      self.rotations,
                                      self.translations,
                                      self.symm_map)

    def adjust_stress(self, atoms, stress):
        stress_3x3 = voigt_6_to_full_3x3_stress(stress)
        stress_3x3[2] = 0.0
        stress[:] = full_3x3_to_voigt_6_stress(stress_3x3)


        raw_stress = voigt_6_to_full_3x3_stress(stress)
        symmetrized_stress = symmetrize_stress(atoms.get_cell(),
                                               atoms.get_reciprocal_cell().T,
                                               raw_stress, self.rotations)
        stress[:] = full_3x3_to_voigt_6_stress(symmetrized_stress)
