from ase.structure import molecule
from ase.calculators.calculator import get_calculator


required = {'abinit': dict(ecut=200, toldfe=0.0001)}


def aep1_test(name, k=True, dft=True, restart=True):
    Calculator = get_calculator(name)
    par = required.get(name, {})
    if restart:
        par['label'] = name
    if dft:
        par['xc'] = 'LDA'
    calc = Calculator(**par)
    h2 = molecule('H2', calculator=calc)
    h2.center(vacuum=2.0)
    e2 = h2.get_potential_energy()
    if dft:
        calc.set(xc='PBE')
        e2pbe = h2.get_potential_energy()
    h1 = h2.copy()
    del h1[1]
    h1.set_initial_magnetic_moments([1])
    h1.calc = calc
    if dft:
        e1pbe = h1.get_potential_energy()
        calc.set(xc='LDA')
    e1 = h1.get_potential_energy()
    try:
        m1 = h1.get_magnetic_moment()
    except NotImplementedError:
        pass
    else:
        print m1
    print(2 * e1 - e2)
    if dft:
        print(2 * e1pbe - e2pbe)
    if restart:
        calc = Calculator(name)
        print calc.parameters, calc.results, calc.state
        assert not calc.calculation_required(h1, ['energy'])
        h1 = calc.get_atoms()
        print h1.get_potential_energy()
        dir = 'dir/' + name + '-h1'
        par['label'] = dir
        calc = Calculator(atoms=h1, **par)
        print h1.get_potential_energy()
        print Calculator.read_atoms(dir).get_potential_energy()

aep1_test('abinit')
