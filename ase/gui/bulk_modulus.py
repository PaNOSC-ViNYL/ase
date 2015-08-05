import numpy as np

from ase.utils.eos import EquationOfState


def bulk_modulus(images):
    v = np.array([abs(np.linalg.det(A)) for A in images.A])
    EquationOfState(v, images.E).plot()
