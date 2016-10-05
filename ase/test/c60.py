import numpy as np

from ase.build import molecule

from ase.utils.ff import Morse, Angle, Dihedral, VdW, rel_pos_pbc
from ase.calculators.ff import ForceField

from ase.optimize.precon.lbfgs import PreconLBFGS
from ase.optimize.precon import FF

a = molecule('C60')
a.set_cell(50.0*np.identity(3))

# force field parameters for fulleren, Z. Berkai at al. Energy Procedia, 74, 2015, 59-64
cutoff = 1.5
morse_D = 6.1322
morse_alpha = 1.8502
morse_r0 = 1.4322
angle_k = 10.0
angle_a0 = np.deg2rad(120.0)
dihedral_k = 0.346
vdw_epsilonij = 0.0115
vdw_rminij = 3.4681

neighbor_list = [[] for _ in range(len(a))]
morses = []
angles = []
dihedrals = []
vdws = []

# create neighbor list and list of morse interactions
for i in range(len(a)):
    for j in range(i+1, len(a)):
        rij = rel_pos_pbc(a, i, j)
        dij = np.linalg.norm(rij)
        if dij > cutoff:
            continue
        morses.append(Morse(atomi=i, atomj=j, D=morse_D, alpha=morse_alpha, r0=morse_r0))
        neighbor_list[i].append(j)
        neighbor_list[j].append(i)

# create lists of bending and torsion interactions
for i in range(len(a)):
    for jj in range(len(neighbor_list[i])):
        j = neighbor_list[i][jj]
        for kk in range(jj+1, len(neighbor_list[i])):
            k = neighbor_list[i][kk]
            angles.append(Angle(atomi=j, atomj=i, atomk=k, k=angle_k, a0=angle_a0, cos=True))
            for ll in range(kk+1, len(neighbor_list[i])):
                l = neighbor_list[i][ll]
                dihedrals.append(Dihedral(atomi=j, atomj=i, atomk=k, atoml=l, k=dihedral_k))

# create list of van der Waals interactions
for i in range(len(a)):
    morses_angles = []
    morses_angles.extend(neighbor_list[i])
    for n in range(len(neighbor_list[i])):
        morses_angles.extend(neighbor_list[neighbor_list[i][n]])
    for j in range(i+1, len(a)):
        if j not in morses_angles:
            vdws.append(VdW(atomi=i, atomj=j, epsilonij=vdw_epsilonij, rminij=vdw_rminij))

# set up ForceField calculator
calc = ForceField(morses=morses, angles=angles, dihedrals=dihedrals, vdws=vdws)

a1 = a.copy()
a1.set_calculator(calc)
a1.rattle(0.05)

# geometry optimisation without preconditioner
opt = PreconLBFGS(a1, use_armijo=True, precon='ID')
opt.run(fmax=0.0001)
e1 = a1.get_potential_energy()

a2 = a.copy()
a2.set_calculator(calc)
a2.rattle(0.05)

# geometry optimisation with FF based preconditioner
precon = FF(morses=morses, angles=angles, dihedrals=dihedrals)

opt = PreconLBFGS(a2, use_armijo=True, precon=precon)
opt.run(fmax=0.0001)
e2 = a2.get_potential_energy()

assert e1 - 17.238525 < 1e-6
assert e2 - 17.238525 < 1e-6
