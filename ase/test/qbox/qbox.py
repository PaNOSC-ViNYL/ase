"""Tests related to QBOX"""

import os
import numpy as np

from ase import Atoms
from ase.io import qbox
from ase.io import formats
import ase

test_file = os.path.join(os.path.dirname(ase.__file__), 'test', 'qbox',
                         'test.xml')


def read_output():
    """Test reading the output file"""

    # Read only one frame
    atoms = qbox.read_qbox(test_file)

    assert isinstance(atoms, Atoms)
    assert np.allclose(atoms.cell, np.diag([16, 16, 16]))

    assert len(atoms) == 4
    assert np.allclose(atoms[0].position,
                       [3.70001108, -0.00000000, -0.00000003],
                       atol=1e-7)  # Last frame
    assert np.allclose(atoms.get_velocities()[2],
                       [-0.00000089, -0.00000000, -0.00000000],
                       atol=1e-9)  # Last frame
    assert np.allclose(atoms.get_forces()[3],
                       [-0.00000026, -0.01699708, 0.00000746],
                       atol=1e-7)  # Last frame
    assert np.isclose(-15.37294664, atoms.get_potential_energy())
    assert np.allclose(atoms.get_stress(),
                       [-0.40353661, -1.11698386, -1.39096418,
                        0.00001786, -0.00002405, -0.00000014])

    # Read all the frames
    atoms = qbox.read_qbox(test_file, slice(None))

    assert isinstance(atoms, list)
    assert len(atoms) == 5

    assert len(atoms[1]) == 4
    assert np.allclose(atoms[1][0].position,
                       [3.70001108, -0.00000000, -0.00000003],
                       atol=1e-7)  # 2nd frame
    assert np.allclose(atoms[1].get_forces()[3],
                       [-0.00000029, -0.01705361, 0.00000763],
                       atol=1e-7)  # 2nd frame


def test_format():
    """Make sure the `formats.py` operations work"""

    atoms = formats.read(test_file)
    assert len(atoms) == 4

    atoms = formats.read(test_file, index=slice(None), format='qbox')
    assert len(atoms) == 5


if __name__ == '__main__':
    read_output()
    test_format()
