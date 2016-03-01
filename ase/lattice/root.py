from math import atan2, ceil, cos, sin, log10
import numpy as np


def root_surface(primitive_slab, root, cell_vectors=None, swap_alpha=False,
                 return_valid=False, eps=1e-8):
    """This script allows a surface to be maniupulated to repeat with a 
    special root cut cell.  The general use of this is to make cells with 
    a specific cell geometry, but a nonstandard number of repetitions in
    the cell.  Without using a tool like this, it would be impossible to
    trivially make a fcc111 cell with 13 atoms on each layer, while still
    preserving the geometry of the primitive cell.

    *primitive cell* should be a primitive 2d cell of your slab, repeated 
    as needed in the z direction
    *root* should be determined using an analysis tool such as the 
    root_surface_analysis function.
    *cell_vectors* is a manual override for the detected cell
    *swap_alpha* swaps the alpha angle of the cell
    *eps* is a precision value as this relies on floating point precision, 
    adjust when more accurate cells are needed but the default seems to 
    work with all tested cells"""

    atoms = primitive_slab.copy()
    # If cell_vectors is not given, try to guess from the atoms
    # Normalize the x axis to a distance of 1, and use the cell
    # We ignore the z axis because this code cannot handle it
    if cell_vectors is None:
        xscale = np.linalg.norm(atoms._cell[0][0:2])
        xx, xy = atoms._cell[0][0:2] / xscale
        yx, yy = atoms._cell[1][0:2] / xscale
        cell_vectors = [[xx, xy], [yx, yy]]

    # Make (0, 0) corner's angle flip from acute to obtuse or
    # obtuse to acute with a small trick
    if swap_alpha:
        cell_vectors[1][0] *= -1

    # Manipulate the cell vectors to find the best search zone and
    # cast to numpy array.
    cell_vectors = np.array(cell_vectors)
    cell_vectors_mag = map(np.linalg.norm, cell_vectors)
    cell_search = map(lambda x: int(ceil(float(root) / float(x))),
                      cell_vectors_mag)

    # Make these variables in function scope
    # x,  y  = Raw grid point
    # tx, ty = Transformed grid point
    x, y, tx, ty = 0, 0, 0, 0

    # Returns valid roots that are found in the given search
    # space.  To find more, use a higher root.
    if return_valid:
        valid = set()
        for x in range(cell_search[0]):
            for y in range(cell_search[1]):
                if x == y == 0:
                    continue
                vect = (cell_vectors[0] * x) + (cell_vectors[1] * y)
                dist = (vect ** 2).sum()
                valid.add(dist)
        return sorted(list(valid))

    # Calculate square distances and break when appropriate
    for x in range(cell_search[0]):
        for y in range(cell_search[1]):
            if x == y == 0:
                continue
            vect = (cell_vectors[0] * x) + (cell_vectors[1] * y)
            dist = (vect ** 2).sum()
            if abs(dist - root) <= eps:
                tx, ty = vect
                break
        else:
            continue
        break
    else:
        # A root cell could not be found for this combination
        raise RuntimeError("Can't find a root cell of {0} in [{1}, {2}]".
                           format(root, cell_vectors[0], cell_vectors[1]))

    tmag = np.linalg.norm((tx, ty))
    root_angle = atan2(ty, tx)
    cell_scale = tmag / cell_vectors_mag[0]
    root_rotation = [[cos(root_angle), -sin(root_angle)],
                     [sin(root_angle), cos(root_angle)]]
    cell = map(lambda x: np.dot(x, root_rotation) * cell_scale, cell_vectors)

    def remove_doubles(atoms, shift=True):
        if shift:
            atoms.translate((0.0007, 0.0008, 0.0009))
        atoms.set_scaled_positions(atoms.get_scaled_positions())
        valid = [0]
        for x in range(len(atoms)):
            for ypos, y in enumerate(valid):
                xa = atoms[x].position
                ya = atoms[y].position
                if np.linalg.norm(xa - ya) < eps:
                    break
            else:
                valid.append(x)
        del atoms[[i for i in range(len(atoms)) if i not in valid]]
        if shift:
            atoms.translate((-0.0007, -0.0008, -0.0009))

    atoms_cell_mag = map(np.linalg.norm, np.array(atoms._cell[0:2, 0:2]))
    cell_vect_mag = map(np.linalg.norm, np.array(cell_vectors))
    cell_scale = np.divide(atoms_cell_mag, cell_vect_mag)
    atoms *= (cell_search[0], cell_search[1], 1)
    atoms._cell[0:2, 0:2] = cell * cell_scale
    remove_doubles(atoms, shift=False)
    remove_doubles(atoms, shift=True)

    def rot(vector, angle):
        return [(vector[0] * cos(angle)) - (vector[1] * sin(angle)),
                (vector[0] * sin(angle)) + (vector[1] * cos(angle))]
    angle = -atan2(atoms._cell[0][1], atoms._cell[0][0])
    atoms._cell[0][0:2] = rot(atoms._cell[0][0:2], angle)
    atoms._cell[1][0:2] = rot(atoms._cell[1][0:2], angle)
    for atom in atoms:
        atom.position[0:2] = rot(atom.position[0:2], angle)
    atoms.center()

    atoms.positions = np.around(atoms.positions, decimals=int(-log10(eps)))
    ind = np.lexsort(
        (atoms.positions[:, 0], atoms.positions[:, 1], atoms.positions[:, 2],))
    return atoms[ind]

def root_surface_analysis(primitive_slab, root, cell_vectors=None):
    """This is a tool to analyze a slab and look for valid roots that exist,
       without using this, nontrivial cells may be difficult to find."""
    atoms = primitive_slab
    # If cell_vectors is not given, try to guess from the atoms
    # Normalize the x axis to a distance of 1, and use the cell
    # We ignore the z axis because this code cannot handle it
    if cell_vectors is None:
        xscale = np.linalg.norm(atoms._cell[0][0:2])
        xx, xy = atoms._cell[0][0:2] / xscale
        yx, yy = atoms._cell[1][0:2] / xscale
        cell_vectors = [[xx, xy], [yx, yy]]

    # Manipulate the cell vectors to find the best search zone and
    # cast to numpy array.
    cell_vectors = np.array(cell_vectors)
    cell_vectors_mag = map(np.linalg.norm, cell_vectors)
    cell_search = map(lambda x: int(ceil(float(root) / float(x))),
                      cell_vectors_mag)

    # Returns valid roots that are found in the given search
    # space.  To find more, use a higher root.
    if return_valid:
        valid = set()
        for x in range(cell_search[0]):
            for y in range(cell_search[1]):
                if x == y == 0:
                    continue
                vect = (cell_vectors[0] * x) + (cell_vectors[1] * y)
                dist = (vect ** 2).sum()
                valid.add(dist)
        return sorted(list(valid))

