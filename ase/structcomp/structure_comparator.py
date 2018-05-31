from __future__ import print_function
import copy
import numpy as np
from scipy.spatial import cKDTree as KDTree
from ase import Atoms
from ase.build import tools as asetools

try:
    import pystructcomp_cpp as pycpp
    # The C++ version is available at:
    # https://github.com/davidkleiven/StructureCompare
    has_cpp_version = True
except:
    has_cpp_version = False


class SymmetryEquivalenceCheck(object):
    """Compare two structures to determine if they are symmetry equivalent.
    
    Based on the recipe from Comput. Phys. Commun. 183, 690–697 (2012).

    Attributes:
    ==========
    angle_tol: float
        angle tolerance for the lattice vectors in degrees

    ltol: float
        relative tolerance for the length of the lattice vectors (per atom)

    stol: float
        position tolerance for the site comparison in units of
        (V/N)^(1/3) (average length between atoms)

    scale_volume: bool
        if True the volumes of the two structures are scaled to be equal
    """
    def __init__(self, angle_tol=1.0, ltol=0.05, stol=0.05,
                 scale_volume=False):
        self.s1 = None
        self.s2 = None
        self.angle_tol = angle_tol * np.pi / 180.0  # convert to radians
        self.scale_volume = scale_volume
        self.stol = stol
        self.ltol = ltol
        self.position_tolerance = 0.0
        self.use_cpp_version = has_cpp_version

    def _niggli_reduce(self):
        """Reduce to niggli cells.

        Reduce the two atoms to niggli cells, then rotates the niggli cells to
        the so called "standard" orientation with one lattice vector along the
        x-axis.
        """
        asetools.niggli_reduce(self.s1)
        asetools.niggli_reduce(self.s2)
        self._standarize_cell(self.s1)
        self._standarize_cell(self.s2)

    def _standarize_cell(self, atoms):
        """Rotate the first vector such that it points along the x-axis."""
        cell = atoms.get_cell().T
        total_rot_mat = np.eye(3)
        v1 = cell[:, 0]
        l1 = np.sqrt(v1[0]**2 + v1[2]**2)
        angle = np.abs(np.arcsin(v1[2] / l1))
        if (v1[0] < 0.0 and v1[2] > 0.0):
            angle = np.pi - angle
        elif (v1[0] < 0.0 and v1[2] < 0.0):
            angle = np.pi + angle
        elif (v1[0] > 0.0 and v1[2] < 0.0):
            angle = -angle
        ca = np.cos(angle)
        sa = np.sin(angle)
        rotmat = np.array([[ca, 0.0, sa], [0.0, 1.0, 0.0], [-sa, 0.0, ca]])
        total_rot_mat = rotmat.dot(total_rot_mat)
        cell = rotmat.dot(cell)

        v1 = cell[:, 0]
        l1 = np.sqrt(v1[0]**2 + v1[1]**2)
        angle = np.abs(np.arcsin(v1[1] / l1))
        if (v1[0] < 0.0 and v1[1] > 0.0):
            angle = np.pi - angle
        elif (v1[0] < 0.0 and v1[1] < 0.0):
            angle = np.pi + angle
        elif (v1[0] > 0.0 and v1[1] < 0.0):
            angle = -angle
        ca = np.cos(angle)
        sa = np.sin(angle)
        rotmat = np.array([[ca, sa, 0.0], [-sa, ca, 0.0], [0.0, 0.0, 1.0]])
        total_rot_mat = rotmat.dot(total_rot_mat)
        cell = rotmat.dot(cell)

        # Rotate around x axis such that the second vector is in the xy plane
        v2 = cell[:, 1]
        l2 = np.sqrt(v2[1]**2 + v2[2]**2)
        angle = np.abs(np.arcsin(v2[2] / l2))
        if (v2[1] < 0.0 and v2[2] > 0.0):
            angle = np.pi - angle
        elif (v2[1] < 0.0 and v2[2] < 0.0):
            angle = np.pi + angle
        elif (v2[1] > 0.0 and v2[2] < 0.0):
            angle = -angle
        ca = np.cos(angle)
        sa = np.sin(angle)
        rotmat = np.array([[1.0, 0.0, 0.0], [0.0, ca, sa], [0.0, -sa, ca]])
        total_rot_mat = rotmat.dot(total_rot_mat)
        cell = rotmat.dot(cell)

        atoms.set_cell(cell.T)
        atoms.set_positions(total_rot_mat.dot(atoms.get_positions().T).T)
        atoms.wrap(pbc=[1, 1, 1])
        return atoms

    def _get_element_count(self):
        """Count the number of elements in each of the structures."""
        elem1 = {}
        elem2 = {}

        for atom in self.s1:
            if atom.symbol in elem1.keys():
                elem1[atom.symbol] += 1
            else:
                elem1[atom.symbol] = 1

        for atom in self.s2:
            if atom.symbol in elem2.keys():
                elem2[atom.symbol] += 1
            else:
                elem2[atom.symbol] = 1
        return elem1, elem2

    def _has_same_elements(self):
        """Check if two structures have same elements."""
        elem1, elem2 = self._get_element_count()
        return elem1 == elem2

    def _has_same_angles(self):
        """Check that the Niggli unit vectors has the same internal angles."""
        cell1 = self.s1.get_cell().T
        cell2 = self.s2.get_cell().T

        # Normalize each vector
        for i in range(3):
            cell1[:, i] /= np.sqrt(np.sum(cell1[:, i]**2))
            cell2[:, i] /= np.sqrt(np.sum(cell2[:, i]**2))
        dot1 = cell1.T.dot(cell1)
        dot2 = cell2.T.dot(cell2)

        # Extract only the relevant dot products
        dot1 = [dot1[0, 1], dot1[0, 2], dot1[1, 2]]
        dot2 = [dot2[0, 1], dot2[0, 2], dot2[1, 2]]

        # Convert to angles
        ang1 = [np.arccos(scalar_prod) * 180.0 / np.pi for scalar_prod in dot1]
        ang2 = [np.arccos(scalar_prod) * 180.0 / np.pi for scalar_prod in dot2]

        for i in range(3):
            closestIndex = np.argmin(np.abs(np.array(ang2) - ang1[i]))
            if np.abs(ang2[closestIndex] - ang1[i]) < self.angle_tol:
                # Remove the entry that matched
                del ang2[closestIndex]
            else:
                return False
        return True

    def _scale_volumes(self):
        """Scale the cell of s1 to have the same volume as s2."""
        v1 = np.linalg.det(self.s1.get_cell())
        v2 = np.linalg.det(self.s2.get_cell())

        # Scale the cells
        cell1 = self.s1.get_cell()
        coordinate_scaling = (v2 / v1)**(1.0 / 3.0)
        cell1 *= coordinate_scaling
        self.s1.set_positions(self.s1.get_positions() * coordinate_scaling)
        self.s1.set_cell(cell1)

    def _has_same_volume(self):
        vol1 = self.s1.get_volume()
        vol2 = self.s2.get_volume()
        return np.abs(vol1 - vol2) < 1E-5

    def compare(self, s1, s2):
        """Compare the two structures.

        Return *True* if the two structures are equivalent, *False* otherwise.

        Arguments:
        =========
        s1, s2: Atoms objects
        """
        self.s1 = s1.copy()
        self.s2 = s2.copy()

        vol = self.s1.get_volume()
        self.position_tolerance = self.stol * (vol / len(self.s2))**(1.0 / 3.0)

        if len(s1) != len(s2):
            return False

        if not self._has_same_elements():
            return False

        self._niggli_reduce()
        if not self._has_same_angles():
            return False

        if self.scale_volume:
            self._scale_volumes()

        if not self._has_same_volume():
            return False

        if self.use_cpp_version:
            return self._compare_cpp()

        matrices, translations = self._get_rotation_reflection_matrices()
        return self._positions_match(matrices, translations)

    def _compare_cpp(self):
        """Compare the two structures using the C++ version."""
        atoms1_ref, atoms2_ref = \
            self._extract_positions_of_least_frequent_element()
        cell = self.s1.get_cell().T

        # Additional vector that is added to make sure that
        # there always is an atom at the origin
        delta_vec = 1E-6 * (cell[:, 0] + cell[:, 1] + cell[:, 2])

        # Put on of the least frequent elements of structure 2 at the origin
        translation = atoms2_ref.get_positions()[0, :] - delta_vec
        atoms2_ref.set_positions(atoms2_ref.get_positions() - translation)
        atoms2_ref.wrap(pbc=[1, 1, 1])
        self.s2.set_positions(self.s2.get_positions() - translation)
        self.s2.wrap(pbc=[1, 1, 1])
        translation = atoms1_ref.get_positions()[0, :] - delta_vec

        sc_atom_search = atoms1_ref * (3, 3, 3)
        sc_pos = atoms1_ref.get_positions()
        # Translate by one cell diagonal
        sc_pos += cell[:, 0] + cell[:, 1] + cell[:, 2]

        sc_pos_search = sc_atom_search.get_positions()

        translation = sc_pos[0, :] - delta_vec

        new_sc_pos = sc_pos_search - translation
        sc_atom_search.set_positions(new_sc_pos)
        sc_atom_search.wrap()
        exp2 = self._expand(self.s2)
        symb1 = [atom.symbol for atom in self.s1]
        symb_exp = [atom.symbol for atom in exp2]
        return pycpp.compare(self, self.s1, exp2, atoms1_ref, sc_atom_search,
                             symb1, symb_exp)

    def _get_least_frequent_element(self):
        """Return the symbol of the least frequent element."""
        elem1, elem2 = self._get_element_count()
        assert elem1 == elem2
        minimum_value = 2 * len(self.s1)  # Set the value to a large value
        least_freq_element = "X"
        for key, value in elem1.items():
            if value < minimum_value:
                least_freq_element = key
                minimum_value = value

        if least_freq_element == "X":
            msg = "Did not manage to find the least frequent element"
            raise ValueError(msg)
        return least_freq_element

    def _extract_positions_of_least_frequent_element(self):
        """Extract a dictionary of positions of each element."""
        least_freq_element = self._get_least_frequent_element()

        atoms1 = None
        atoms2 = None

        for i in range(len(self.s1)):
            symbol = self.s1[i].symbol
            if symbol == least_freq_element:
                pos = self.s1.get_positions(wrap=True)[i, :]
                if atoms1 is None:
                    atoms1 = Atoms(symbol, positions=[pos])
                else:
                    atoms1.extend(Atoms(symbol, positions=[pos]))

            symbol = self.s2[i].symbol
            if symbol == least_freq_element:
                pos = self.s2.get_positions(wrap=True)[i, :]
                if atoms2 is None:
                    atoms2 = Atoms(symbol, positions=[pos])
                else:
                    atoms2.extend(Atoms(symbol, positions=[pos]))
        atoms1.set_cell(self.s1.get_cell())
        atoms2.set_cell(self.s2.get_cell())
        return atoms1, atoms2

    def _positions_match(self, rotation_reflection_matrices, translations):
        """Check if the position and elements match.

        Note that this function changes self.s1 and self.s2 to the rotation and
        translation that matches best. Hence, it is crucial that this function
        is called before the element comparison.
        """
        # Position matching not implemented yet
        pos1_ref = self.s1.get_positions(wrap=True)

        # Expand the reference object
        exp2 = self._expand(self.s2)

        # Build a KD tree to enable fast look-up of nearest neighbours
        tree = KDTree(exp2.get_positions())
        for i in range(translations.shape[0]):
            for matrix in rotation_reflection_matrices:
                pos1 = copy.deepcopy(pos1_ref)
                # Translate
                pos1 -= translations[i, :]

                # Rotate
                pos1 = matrix.dot(pos1.T).T

                # Update the atoms positions
                self.s1.set_positions(pos1)
                self.s1.wrap(pbc=[1, 1, 1])
                if self._elements_match(self.s1, exp2, tree):
                    return True
        return False

    def _expand(self, ref_atoms):
        """Add additional atoms to for atoms close to the cell boundaries.

        This ensures that atoms having crossed the cell boundaries due to
        numerical noise are properly detected.
        """
        expanded_atoms = copy.deepcopy(ref_atoms)

        cell = ref_atoms.get_cell().T
        normal_vectors = [np.cross(cell[:, 0], cell[:, 1]),
                          np.cross(cell[:, 0], cell[:, 2]),
                          np.cross(cell[:, 1], cell[:, 2])]
        normal_vectors = [vec / np.sqrt(np.sum(vec**2)) for vec
                          in normal_vectors]
        positions = ref_atoms.get_positions(wrap=True)
        tol = 0.0001

        for i in range(len(ref_atoms)):
            surface_close = [False, False, False, False, False, False]
            # Face 1
            distance = np.abs(positions[i, :].dot(normal_vectors[0]))
            symbol = ref_atoms[i].symbol
            if distance < tol:
                newpos = positions[i, :] + cell[:, 2]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
                surface_close[0] = True

            # Face 2
            vec = positions[i, :] - cell[:, 2]
            distance = np.abs(vec.dot(normal_vectors[0]))
            if distance < tol:
                newpos = positions[i, :] - cell[:, 2]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
                surface_close[1] = True

            # Face 3
            distance = np.abs(positions[i, :].dot(normal_vectors[1]))
            if distance < tol:
                newpos = positions[i, :] + cell[:, 1]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
                surface_close[2] = True

            # Face 4
            vec = positions[i, :] - cell[:, 1]
            distance = np.abs(vec.dot(normal_vectors[1]))
            if distance < tol:
                newpos = positions[i, :] - cell[:, 1]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
                surface_close[3] = True

            # Face 5
            distance = np.abs(positions[i, :].dot(normal_vectors[2]))
            if distance < tol:
                newpos = positions[i, :] + cell[:, 0]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
                surface_close[4] = True

            # Face 6
            vec = positions[i, :] - cell[:, 0]
            distance = np.abs(vec.dot(normal_vectors[2]))
            if distance < tol:
                newpos = positions[i, :] - cell[:, 0]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
                surface_close[5] = True

            # Take edges into account
            if (surface_close[0] and surface_close[2]):
                newpos = positions[i, :] + cell[:, 1] + cell[:, 2]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
            if (surface_close[0] and surface_close[4]):
                newpos = positions[i, :] + cell[:, 0] + cell[:, 2]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
            if (surface_close[0] and surface_close[3]):
                newpos = positions[i, :] - cell[:, 1] + cell[:, 2]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
            if (surface_close[0] and surface_close[5]):
                newpos = positions[i, :] - cell[:, 0] + cell[:, 2]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
            if (surface_close[3] and surface_close[4]):
                newpos = positions[i, :] + cell[:, 0] - cell[:, 1]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
            if (surface_close[3] and surface_close[1]):
                newpos = positions[i, :] - cell[:, 1] - cell[:, 2]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
            if (surface_close[3] and surface_close[5]):
                newpos = positions[i, :] - cell[:, 1] - cell[:, 0]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
            if (surface_close[1] and surface_close[5]):
                newpos = positions[i, :] - cell[:, 0] - cell[:, 2]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
            if (surface_close[2] and surface_close[5]):
                newpos = positions[i, :] - cell[:, 0] + cell[:, 1]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
            if (surface_close[2] and surface_close[4]):
                newpos = positions[i, :] + cell[:, 0] + cell[:, 1]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
            if (surface_close[1] and surface_close[2]):
                newpos = positions[i, :] + cell[:, 1] - cell[:, 2]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
            if (surface_close[1] and surface_close[4]):
                newpos = positions[i, :] + cell[:, 0] - cell[:, 2]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)

            # Take corners into account
            if (surface_close[0] and surface_close[2] and surface_close[4]):
                newpos = positions[i, :] + cell[:, 0] + cell[:, 1] + cell[:, 2]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
            elif (surface_close[0] and surface_close[3] and surface_close[4]):
                newpos = positions[i, :] + cell[:, 0] - cell[:, 1] + cell[:, 2]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
            elif (surface_close[0] and surface_close[2] and surface_close[5]):
                newpos = positions[i, :] - cell[:, 0] + cell[:, 1] + cell[:, 2]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
            elif (surface_close[0] and surface_close[3] and surface_close[5]):
                newpos = positions[i, :] - cell[:, 0] - cell[:, 1] + cell[:, 2]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
            elif (surface_close[1] and surface_close[3] and surface_close[4]):
                newpos = positions[i, :] + cell[:, 0] - cell[:, 1] - cell[:, 2]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
            elif (surface_close[1] and surface_close[3] and surface_close[5]):
                newpos = positions[i, :] - cell[:, 0] - cell[:, 1] - cell[:, 2]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
            elif (surface_close[1] and surface_close[2] and surface_close[4]):
                newpos = positions[i, :] + cell[:, 0] + cell[:, 1] - cell[:, 2]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)
            elif (surface_close[1] and surface_close[2] and surface_close[5]):
                newpos = positions[i, :] - cell[:, 0] + cell[:, 1] - cell[:, 2]
                newAtom = Atoms(symbol, positions=[newpos])
                expanded_atoms.extend(newAtom)

        return expanded_atoms

    def _elements_match(self, s1, s2, kdtree):
        """Check if all the elements in s1 match the corresponding position in s2

        NOTE: The unit cells may be in different octants
        Hence, try all cyclic permutations of x,y and z
        """
        pos1 = s1.get_positions()
        for order in range(1):
            all_match = True
            for i in range(len(s1)):
                s1pos = np.zeros(3)
                s1pos[0] = pos1[i, order]
                s1pos[1] = pos1[i, (order + 1) % 3]
                s1pos[2] = pos1[i, (order + 2) % 3]
                dist, closest = kdtree.query(s1pos)
                if (s1[i].symbol != s2[closest].symbol or
                        dist > self.position_tolerance):
                    all_match = False
                    break
            if all_match:
                return True
        return False

    def _get_rotation_reflection_matrices(self):
        """Compute candidates for the transformation matrix."""
        atoms1_ref, atoms2_ref = \
            self._extract_positions_of_least_frequent_element()
        cell = self.s1.get_cell().T
        angle_tol = self.angle_tol

        # Additional vector that is added to make sure that
        # there always is an atom at the origin
        delta_vec = 1E-6 * (cell[:, 0] + cell[:, 1] + cell[:, 2])

        # Put on of the least frequent elements of structure 2 at the origin
        translation = atoms2_ref.get_positions()[0, :] - delta_vec
        atoms2_ref.set_positions(atoms2_ref.get_positions() - translation)
        atoms2_ref.wrap(pbc=[1, 1, 1])
        self.s2.set_positions(self.s2.get_positions() - translation)
        self.s2.wrap(pbc=[1, 1, 1])

        # Store three reference vectors
        ref_vec = atoms2_ref.get_cell().T
        ref_vec_lengths = np.sqrt(np.sum(ref_vec**2, axis=0))

        canditate_trans_mat = []

        # Compute ref vec angles
        angle12_ref = np.arccos(ref_vec[:, 0].dot(ref_vec[:, 1]) /
                                (ref_vec_lengths[0] * ref_vec_lengths[1]))
        if angle12_ref > np.pi / 2.0:
            angle12_ref = np.pi - angle12_ref
        angle13_ref = np.arccos(ref_vec[:, 0].dot(ref_vec[:, 2]) /
                                (ref_vec_lengths[0] * ref_vec_lengths[2]))
        if angle13_ref > np.pi / 2.0:
            angle13_ref = np.pi - angle13_ref
        angle23_ref = np.arccos(ref_vec[:, 1].dot(ref_vec[:, 2]) /
                                (ref_vec_lengths[1] * ref_vec_lengths[2]))
        if angle23_ref > np.pi / 2.0:
            angle23_ref = np.pi - angle23_ref

        sc_atom_search = atoms1_ref * (3, 3, 3)
        sc_pos = atoms1_ref.get_positions()
        # Translate by one cell diagonal
        sc_pos += cell[:, 0] + cell[:, 1] + cell[:, 2]
        sc_pos_search = sc_atom_search.get_positions()

        candidate_vecs = [[], [], []]
        translation = sc_pos[0, :] - delta_vec

        new_sc_pos = sc_pos_search-translation
        lengths = np.sqrt(np.sum(new_sc_pos**2, axis=1))
        for l in range(1, len(lengths)):
            for k in range(3):
                if (np.abs(lengths[l] - ref_vec_lengths[k]) <
                        self.ltol * lengths[l] / len(self.s1)):
                    candidate_vecs[k].append(new_sc_pos[l, :])

        # Check angles
        refined_candidate_list = [[], [], []]

        for v1 in candidate_vecs[0]:
            for v2 in candidate_vecs[1]:
                if np.allclose(v1, v2, atol=1E-3):
                    continue
                v1len = np.sqrt(np.sum(v1**2))
                v2len = np.sqrt(np.sum(v2**2))
                angle12 = np.arccos(v1.dot(v2) / (v1len * v2len))
                if angle12 > np.pi / 2.0:
                    angle12 = np.pi - angle12
                for v3 in candidate_vecs[2]:
                    if (np.allclose(v1, v3, atol=1E-3) or
                            np.allclose(v2, v3, atol=1E-3)):
                        continue
                    v3len = np.sqrt(np.sum(v3**2))
                    angle13 = np.arccos(v1.dot(v3) / (v1len * v3len))
                    if angle13 > np.pi / 2.0:
                        angle13 = np.pi-angle13

                    angle23 = np.arccos(v2.dot(v3) / (v2len * v3len))
                    if angle23 > np.pi / 2.0:
                        angle23 = np.pi - angle23

                    if (np.abs(angle12 - angle12_ref) < angle_tol and
                            np.abs(angle13 - angle13_ref) < angle_tol and
                            np.abs(angle23 - angle23_ref) < angle_tol):
                        refined_candidate_list[0].append(v1)
                        refined_candidate_list[1].append(v2)
                        refined_candidate_list[2].append(v3)

        # Compute rotation/reflection
        for v1, v2, v3 in zip(refined_candidate_list[0],
                              refined_candidate_list[1],
                              refined_candidate_list[2]):
            T = np.zeros((3, 3))
            T[:, 0] = v1
            T[:, 1] = v2
            T[:, 2] = v3
            R = ref_vec.dot(np.linalg.inv(T))
            canditate_trans_mat.append(R)
        return canditate_trans_mat, atoms1_ref.get_positions()
