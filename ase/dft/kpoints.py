import numpy as np


def monkhorst_pack(size):
    """Construct a uniform sampling of k-space of given size."""
    if np.less_equal(size, 0).any():
        raise ValueError('Illegal size: %s' % list(size))
    kpts = np.indices(size).transpose((1, 2, 3, 0)).reshape((-1, 3))
    return (kpts + 0.5) / size - 0.5


def get_monkhorst_shape(kpts, tol=1e-5):
    """Return the number of k-points along each axis of input Monkhorst pack.

    The set of k-points must not have been symmetry reduced.
    """
    nkpts = len(kpts)
    if nkpts == 1:
        return np.ones(3, int)
    
    Nk_c = np.zeros(3, int)
    for c in range(3):
        # Determine increment between kpoints along current axis
        DeltaK = max(np.diff(np.sort(kpts[:, c])))

        # Determine number of kpoints as inverse of distance between kpoints
        if DeltaK > tol:
            Nk_c[c] = int(round(1. / DeltaK))
        else:
            Nk_c[c] = 1
    return Nk_c


def kpoint_convert(cell_cv, skpts_kc=None, ckpts_kv=None):
    """Convert k-points between scaled and cartesian coordinates.

    Given the atomic unit cell, and either the scaled or cartesian k-point
    coordinates, the other is determined.

    The k-point arrays can be either a single point, or a list of points,
    i.e. the dimension k can be empty or multidimensional.
    """
    if ckpts_kv is None:
        icell_cv = 2 * np.pi * np.linalg.inv(cell_cv).T
        return np.dot(skpts_kc, icell_cv)
    elif skpts_kc is None:
        return np.dot(ckpts_kv, cell_cv.T) / (2 * np.pi)
    else:
        raise KeyError('Either scaled or cartesian coordinates must be given.')


def get_bandpath(points, cell, npoints=50):
    """Make a list of kpoints defining the path between the given points.

    points is a list of special IBZ point pairs, e.g.
    >>> points = [L, Gamma, Gamma, X, X, U, K, Gamma]
    These should be given in scaled coordinates.

    cell is the unitcell if the atoms.

    npoints is the approximate desired length of the output kpts list

    The output list point_indices, gives the indices in kpts of the special
    points.
    """
    assert len(points) % 2 == 0
    points = np.asarray(points)
    dists = points[1::2] - points[::2]
    lengths = [np.linalg.norm(d) for d in kpoint_convert(cell, skpts_kc=dists)]
    length = sum(lengths)
    kpts = []
    point_indices = [0]
    for P, d, L in zip(points[::2], dists, lengths):
        for t in np.linspace(0, 1, int(round(L / length * npoints)),
                             endpoint=False):
            kpts.append(P + t * d)
        point_indices.append(len(kpts))
    return np.array(kpts), point_indices


# The following is a list of the critical points in the 1. Brillouin zone
# for some typical crystal structures.
# (In units of the reciprocal basis vectors)
ibz_points = {'cubic': {'Gamma': [0/1., 0/1., 0/1.],
                        'X':     [0/1., 0/2., 1/2.],
                        'R':     [1/2., 1/2., 1/2.],
                        'M':     [0/2., 1/2., 1/2.]},
              'fcc':   {'Gamma': [0/1., 0/1., 0/1.],
                        'X':     [0/1., 1/2., 1/2.],
                        'W':     [1/2., 3/4., 1/4.],
                        'K':     [3/8., 3/8., 3/4.],
                        'U':     [1/4., 5/8., 5/8.],
                        'L':     [1/2., 1/2., 1/2.]},
              'bcc':   {'Gamma': [0/1., 0/1., 0/1.],
                        'H': None,
                        'N': None,
                        'P': None},
              }
