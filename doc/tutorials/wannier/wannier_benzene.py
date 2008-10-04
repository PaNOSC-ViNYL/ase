from gpaw import restart
from ase.dft import Wannier

atoms, calc = restart('benzene.gpw', txt=None)

wan = Wannier(nwannier=18, calc=calc, fixedstates=15, spin=0)
wan.initialize(calc)
wan.localize()
for i in range(wan.nwannier):
    wan.write_cube(calc, i, 'benzene%i.cube' % i)
