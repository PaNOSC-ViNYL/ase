from ase.geometry.cell import *
from ase.build import bulk

def check_single(name, cell):
    c = Cell(cell)
    name1 = get_structure_name(c)
    ok = name.split('@')[0] == name1
    print(name, '-->', name1, 'OK' if ok else 'ERR', c.cellpar())
    assert ok

def check(name, cell):
    cell = Cell(cell).cell
    # Check all three positive permutations:
    check_single(name + '@012', cell[[0, 1, 2]])
    check_single(name + '@201', cell[[2, 0, 1]])
    check_single(name + '@120', cell[[1, 2, 0]])

check('cub', cub(3.3))
check('fcc', fcc(3.4))
check('fcc', bulk('Au').cell)
check('bcc', bcc(3.5))
check('bcc', bulk('Fe').cell)
check('tet', tet(4., 5.))
check('tet', np.diag([4., 5., 5.]))
check('tet', np.diag([5., 4., 5.]))
check('tet', np.diag([5., 5., 4.]))
check('bct', bct(3., 4.))
check('orc', orc(3., 4., 5.))
check('orc', orc(4., 5., 3.))
check('orcf', orcf(4., 5., 3.))
check('orci', orci(5., 6., 2.))
check('orcc', orcc(4., 3., 5.))
check('hex', hex(5., 6.))
check('rhl', rhl(4., 54.))
check('mcl', mcl(2., 3., 4., 62.))
check('mclc', mclc(3., 5., 4., 70.))
check('tri', tri(7., 6., 5., 34., 55., 80.))
rng = np.random.RandomState(17)
for i in range(3):
    check('tri', rng.rand(3, 3))
