from ase.structcomp import SymmetryEquivalenceCheck
from ase.structcomp.structure_comparator import SpgLibNotFoundError
import os
from ase.test import NotAvailable
from ase.build import bulk
from ase import Atoms
from ase.spacegroup import spacegroup, get_spacegroup, crystal
from random import randint
import numpy as np


def get_atoms_with_mixed_elements(crystalstructure="fcc"):
    atoms = bulk("Al", crystalstructure=crystalstructure, a=3.2)
    atoms = atoms * (2, 2, 2)
    symbs = ["Al", "Cu", "Zn"]
    symbols = [symbs[randint(0, len(symbs) - 1)] for _ in range(len(atoms))]
    for i in range(len(atoms)):
        atoms[i].symbol = symbols[i]
    return atoms


def test_compare(comparator):
    s1 = bulk("Al")
    s1 = s1 * (2, 2, 2)
    s2 = bulk("Al")
    s2 = s2 * (2, 2, 2)
    assert comparator.compare(s1, s2)


def test_fcc_bcc(comparator):
    s1 = bulk("Al", crystalstructure="fcc")
    s2 = bulk("Al", crystalstructure="bcc", a=4.05)
    s1 = s1 * (2, 2, 2)
    s2 = s2 * (2, 2, 2)
    assert not comparator.compare(s1, s2)


def test_single_impurity(comparator):
    s1 = bulk("Al")
    s1 = s1 * (2, 2, 2)
    s1[0].symbol = "Mg"
    s2 = bulk("Al")
    s2 = s2 * (2, 2, 2)
    s2[3].symbol = "Mg"
    assert comparator.compare(s1, s2)


def test_translations(comparator):
    s1 = get_atoms_with_mixed_elements()
    s2 = s1.copy()

    xmax = 2.0 * np.max(s1.get_cell().T)
    N = 1
    dx = xmax / N
    pos_ref = s2.get_positions()
    number_of_correctly_identified = 0
    for i in range(N):
        for j in range(N):
            for k in range(N):
                displacement = np.array([dx * i, dx * j, dx * k])
                new_pos = pos_ref + displacement
                s2.set_positions(new_pos)
                if (comparator.compare(s1, s2)):
                    number_of_correctly_identified += 1

    assert number_of_correctly_identified == N**3


def test_rot_60_deg(comparator):
    s1 = get_atoms_with_mixed_elements()
    s2 = s1.copy()
    ca = np.cos(np.pi / 3.0)
    sa = np.sin(np.pi / 3.0)
    matrix = np.array([[ca, sa, 0.0], [-sa, ca, 0.0], [0.0, 0.0, 1.0]])
    s2.set_positions(matrix.dot(s2.get_positions().T).T)
    s2.set_cell(matrix.dot(s2.get_cell().T).T)
    assert comparator.compare(s1, s2)


def test_rot_120_deg(comparator):
    s1 = get_atoms_with_mixed_elements()
    s2 = s1.copy()
    ca = np.cos(2.0 * np.pi / 3.0)
    sa = np.sin(2.0 * np.pi / 3.0)
    matrix = np.array([[ca, sa, 0.0], [-sa, ca, 0.0], [0.0, 0.0, 1.0]])
    s2.set_positions(matrix.dot(s2.get_positions().T).T)
    s2.set_cell(matrix.dot(s2.get_cell().T).T)
    assert comparator.compare(s1, s2)


def test_rotations_to_standard(comparator):
    s1 = Atoms("Al")
    tol = 1E-6
    for i in range(20):
        cell = np.random.rand(3, 3) * 4.0 - 4.0
        s1.set_cell(cell)
        new_cell = comparator._standarize_cell(s1).get_cell().T
        assert abs(new_cell[1, 0]) < tol
        assert abs(new_cell[2, 0]) < tol
        assert abs(new_cell[2, 1]) < tol


def test_point_inversion(comparator):
    s1 = get_atoms_with_mixed_elements()
    s2 = s1.copy()
    s2.set_positions(-s2.get_positions())
    assert comparator.compare(s1, s2)


def test_mirror_plane(comparator):
    s1 = get_atoms_with_mixed_elements(crystalstructure="hcp")
    s2 = s1.copy()
    mat = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]])
    s2.set_positions(mat.dot(s2.get_positions().T).T)
    assert comparator.compare(s1, s2)

    mat = np.array([[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    s2.set_positions(mat.dot(s1.get_positions().T).T)
    assert comparator.compare(s1, s2)

    mat = np.array([[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]])
    s2.set_positions(mat.dot(s1.get_positions().T).T)
    assert comparator.compare(s1, s2)


def test_hcp_symmetry_ops(comparator):
    s1 = get_atoms_with_mixed_elements(crystalstructure="hcp")
    s2 = s1.copy()
    sg = spacegroup.Spacegroup(194)
    cell = s2.get_cell().T
    inv_cell = np.linalg.inv(cell)
    for op in sg.get_rotations():
        s1 = get_atoms_with_mixed_elements(crystalstructure="hcp")
        s2 = s1.copy()
        transformed_op = cell.dot(op).dot(inv_cell)
        s2.set_positions(transformed_op.dot(s1.get_positions().T).T)
        assert comparator.compare(s1, s2)


def test_fcc_symmetry_ops(comparator):
    s1 = get_atoms_with_mixed_elements()
    s2 = s1.copy()
    sg = spacegroup.Spacegroup(225)
    cell = s2.get_cell().T
    inv_cell = np.linalg.inv(cell)
    for op in sg.get_rotations():
        s1 = get_atoms_with_mixed_elements()
        s2 = s1.copy()
        transformed_op = cell.dot(op).dot(inv_cell)
        s2.set_positions(transformed_op.dot(s1.get_positions().T).T)
        assert comparator.compare(s1, s2)


def test_bcc_symmetry_ops(comparator):
    s1 = get_atoms_with_mixed_elements(crystalstructure="bcc")
    s2 = s1.copy()
    sg = spacegroup.Spacegroup(229)
    cell = s2.get_cell().T
    inv_cell = np.linalg.inv(cell)
    for op in sg.get_rotations():
        s1 = get_atoms_with_mixed_elements(crystalstructure="bcc")
        s2 = s1.copy()
        transformed_op = cell.dot(op).dot(inv_cell)
        s2.set_positions(transformed_op.dot(s1.get_positions().T).T)
        assert comparator.compare(s1, s2)


def test_bcc_translation(comparator):
    s1 = get_atoms_with_mixed_elements(crystalstructure="bcc")
    s2 = s1.copy()
    s2.set_positions(s2.get_positions() + np.array([6.0, -2.0, 1.0]))
    assert comparator.compare(s1, s2)


def test_one_atom_out_of_pos(comparator):
    s1 = get_atoms_with_mixed_elements()
    s2 = s1.copy()
    pos = s1.get_positions()
    pos[0, :] += 0.2
    s2.set_positions(pos)
    assert not comparator.compare(s1, s2)


def test_reduce_to_primitive(comparator):
    try:
        # Tell the comparator to reduce to primitive cell
        comparator.to_primitive = True

        atoms1 = crystal(symbols=['V', 'Li', 'O'],
                         basis=[(0.000000, 0.000000, 0.000000),
                                (0.333333, 0.666667, 0.000000),
                                (0.333333, 0.000000, 0.250000)],
                         spacegroup=167,
                         cellpar=[5.123, 5.123, 13.005, 90., 90., 120.],
                         size=[1, 1, 1], primitive_cell=False)

        atoms2 = crystal(symbols=['V', 'Li', 'O'],
                         basis=[(0.000000, 0.000000, 0.000000),
                                (0.333333, 0.666667, 0.000000),
                                (0.333333, 0.000000, 0.250000)],
                         spacegroup=167,
                         cellpar=[5.123, 5.123, 13.005, 90., 90., 120.],
                         size=[1, 1, 1], primitive_cell=True)

        assert comparator.compare(atoms1, atoms2)
    except SpgLibNotFoundError:
        pass

    # Reset the comparator tot its original state
    comparator.to_primitive = False


def run_all_tests(comparator):
    test_compare(comparator)
    test_fcc_bcc(comparator)
    test_single_impurity(comparator)
    test_translations(comparator)
    test_rot_60_deg(comparator)
    test_rot_120_deg(comparator)
    test_rotations_to_standard(comparator)
    test_point_inversion(comparator)
    test_mirror_plane(comparator)
    test_hcp_symmetry_ops(comparator)
    test_fcc_symmetry_ops(comparator)
    test_bcc_symmetry_ops(comparator)
    test_bcc_translation(comparator)
    test_one_atom_out_of_pos(comparator)
    test_reduce_to_primitive(comparator)


comparator = SymmetryEquivalenceCheck()
if comparator.use_cpp_version:
    run_all_tests(comparator)
    comparator.use_cpp_version = False

run_all_tests(comparator)
