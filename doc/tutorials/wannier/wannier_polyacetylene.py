from ase import *
from ase.dft import Wannier
from gpaw import restart

atoms, calc = restart('poly.gpw', txt=None)

# Make wannier functions using (one) extra degree of freedom
wan = Wannier(nwannier=6, calc=calc, fixedenergy=1.0, file='poly.pickle')
wan.localize()
wan.save('poly.pickle')
wan.translate_all_to_cell((0, 0, 0))
## for i in range(wan.nwannier):
##     wan.write_cube(calc, i, 'polyacetylene_%i.cube' % i)

# Print Kohn-Sham and Wannier bandstructures
ef = calc.get_fermi_level()
## fks = open('KSbands.txt', 'w')
## fwan = open('WANbands.txt', 'w')
## for k, kpt_c in enumerate(calc.get_ibz_k_points()):
##     for eps in calc.get_eigenvalues(kpt=k):
##         print >> fks, kpt_c[0], eps - ef
##     for eps in linalg.eigvalsh(wan.get_hamiltonian(calc, k)).real:
##         print >> fwan, kpt_c[0], eps - ef


fnew = open('NEWbands.txt', 'w')
for k in linspace(-.5, .5, 60):
    for eps in linalg.eigvalsh(wan.get_hamiltonian_kpoint([k, 0, 0], calc,
                                                          cutoff=3.0)).real:
        print >> fnew, k, eps - ef
