from ase import *
from ase.dft import Wannier
from gpaw import restart

atoms, calc = restart('poly.gpw', txt=None)

#wan = Wannier(nwannier=5, calc=calc) # Occ only
wan = Wannier(nwannier=6, calc=calc, occupationenergy=1.0) # One EDF
wan.initialize(calc)
wan.localize()
wan.translate_all_to_cell((2, 0, 0))

for i in range(wan.nwannier):
    wan.write_cube(calc, i, 'polyacetylene%i.cube' % i)
