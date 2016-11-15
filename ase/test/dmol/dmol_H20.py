from ase.build import molecule
from ase.calculators.dmol import DMol3, grd_to_cube
from ase.io import read


atoms = molecule('H2O')
calc = DMol3(symmetry='auto',
             spin_polarization='unrestricted',
             charge=0,
             basis='dnp',
             pseudopotential='none',
             functional='pbe',
             scf_density_convergence=1.0e-7,
             plot=['homo','lumo','density'])

atoms.set_calculator(calc)

print '%-25s %12.5f eV' % ('Potential energy', atoms.get_potential_energy())
print
Econtribs = calc.read_energy_contributions()
for key, val in Econtribs.items():
    print '%-25s %12.5f eV' % (key, val)

print 
print 'index | spin  | Eigenvalue | Occupation'
for spin in [0,1]:
    for i, (eig, occ) in enumerate(zip(calc.get_eigenvalues(spin=spin),
                         calc.get_occupations(spin=spin))):
        print '%5d %5d %12.4f %12.4f' % (i, spin, eig, occ)


atoms_dmol = read(calc.label + '.car', format='dmol-car')
grd_to_cube(atoms_dmol, calc.label + '_density.grd', 'H20_density.cube')
grd_to_cube(atoms_dmol, calc.label + '_homo.grd', 'H20_homo.cube')
grd_to_cube(atoms_dmol, calc.label + '_lumo.grd', 'H20_lumo.cube')

# View cube files
# python -m ase.visualize.mlab H20_lumo.cube
