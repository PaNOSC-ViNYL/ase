import numpy as np
from ase.geometry.cell import Cell, bravais
from ase.build import bulk

def check_single(name, cell):
    c = Cell(cell)
    lattice, params = c.bravais()
    name1 = lattice.type
    ok = name.split('@')[0] == name1
    print(name, '-->', name1, 'OK' if ok else 'ERR', c.cellpar())
    assert ok

def check(name, cell):
    cell = Cell(cell).array
    # Check all three positive permutations:
    check_single(name + '@012', cell[[0, 1, 2]])
    check_single(name + '@201', cell[[2, 0, 1]])
    check_single(name + '@120', cell[[1, 2, 0]])

check('cub', bravais['cub'](3.3))
check('fcc', bravais['fcc'](3.4))
check('fcc', bulk('Au').cell)
check('bcc', bravais['bcc'](3.5))
check('bcc', bulk('Fe').cell)
check('tet', bravais['tet'](4., 5.))
check('tet', np.diag([4., 5., 5.]))
check('tet', np.diag([5., 4., 5.]))
check('tet', np.diag([5., 5., 4.]))
check('bct', bravais['bct'](3., 4.))
check('orc', bravais['orc'](3., 4., 5.))
check('orc', bravais['orc'](4., 5., 3.))
check('orcf', bravais['orcf'](4., 5., 3.))
check('orci', bravais['orci'](5., 6., 2.))
check('orcc', bravais['orcc'](4., 3., 5.))
check('hex', bravais['hex'](5., 6.))
check('rhl', bravais['rhl'](4., 54.))
check('mcl', bravais['mcl'](2., 3., 4., 62.))
check('mclc', bravais['mclc'](3., 5., 4., 70.))
check('tri', bravais['tri'](7., 6., 5., 34., 55., 80.))
rng = np.random.RandomState(17)
for i in range(3):
    tmpcell = rng.rand(3, 3)
    tmpcell *= np.sign(np.linalg.det(tmpcell))
    check('tri', tmpcell)
