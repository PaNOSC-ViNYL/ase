from ase import *

a = 2.70
c = 1.59 * a
h = 1.85
d = 1.10

# XXX Conversion of an old ASE script using Ru.  Geometric data still
# correspond to Ru.
slab = Atoms('2Cu', [(0., 0., 0.), (1/3., 1/3., -.5*c)], tags=(0, 1),
             pbc=(1, 1, 0))

slab.set_cell([(a, 0, 0),
               (a / 2, 3**0.5 * a / 2, 0),
               (0, 0, 1)])
slab = slab.repeat((4, 4, 1))
slab.set_calculator(EMT())
e_slab = slab.get_potential_energy()

molecule = Atoms('2N', positions=[(0., 0., h), (0., 0., h + d)])
molecule.set_calculator(EMT())
e_N2 = molecule.get_potential_energy()

slab.extend(molecule)
constraint = FixAtoms(mask=[a.symbol != 'N' for a in slab])
slab.set_constraint(constraint)
dyn = QuasiNewton(slab)
dyn.run(fmax=0.05)

print 'Adsorption energy:', e_slab + e_N2 - slab.get_potential_energy()
