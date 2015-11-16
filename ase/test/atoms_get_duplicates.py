from ase import Atoms

at = Atoms('H5', positions=[[0., 0., 0.],
                             [1., 0., 0.],
                             [1.01, 0, 0],
                             [3, 2.2, 5.2],
                             [0.1, -0.01, 0.1]])

dups = at.get_duplicate_atoms()
assert all((dups == [[1, 2]]).tolist()) is True

dups = at.get_duplicate_atoms(cutoff=0.2)
assert all((dups == [[0, 4], [1, 2]]).tolist()) is True

at.get_duplicate_atoms(delete=True)
assert len(at) == 4
