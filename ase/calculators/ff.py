from __future__ import division

import numpy as np

from ase.calculators.calculator import Calculator
from ase.utils import ff


class ForceField(Calculator):
    implemented_properties = ['energy', 'forces', 'hessian']
    nolabel = True

    def __init__(self, morses=None, bonds=None, angles=None, dihedrals=None, vdws=None, coulombs=None, **kwargs):
        Calculator.__init__(self, **kwargs)
        if morses is None and bonds is None and angles is None and dihedrals is None and vdws is None and coulombs is None:
           raise ImportError('At least one of morses, bonds, angles, dihedrals, vdws or coulombs lists must be defined!')
        if morses is None:
            self.morses = []
        else:
            self.morses = morses
        if bonds is None:
            self.bonds = []
        else:
            self.bonds = bonds
        if angles is None:
            self.angles = []
        else:
            self.angles = angles
        if dihedrals is None:
            self.dihedrals = []
        else:
            self.dihedrals = dihedrals
        if vdws is None:
            self.vdws = []
        else:
            self.vdws = vdws
        if coulombs is None:
            self.coulombs = []
        else:
            self.coulombs = coulombs

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        if system_changes:
            if 'energy' in self.results:
                del self.results['energy']
            if 'forces' in self.results:
                del self.results['forces']
            if 'hessian' in self.results:
                del self.results['hessian']
        if 'energy' not in self.results:
            energy = 0.0
            for morse in self.morses:
                i, j, e = ff.get_morse_potential_value(atoms, morse)
                energy += e
            for bond in self.bonds:
                i, j, e = ff.get_bond_potential_value(atoms, bond)
                energy += e
            for angle in self.angles:
                i, j, k, e = ff.get_angle_potential_value(atoms, angle)
                energy += e
            for dihedral in self.dihedrals:
                i, j, k, l, e = ff.get_dihedral_potential_value(atoms, dihedral)
                energy += e
            for vdw in self.vdws:
                i, j, e = ff.get_vdw_potential_value(atoms, vdw)
                energy += e
            for coulomb in self.coulombs:
                i, j, e = ff.get_coulomb_potential_value(atoms, coulomb)
                energy += e
            self.results['energy'] = energy
        if 'forces' not in self.results:
            forces = np.zeros(3*len(atoms))
            for morse in self.morses:
                i, j, g = ff.get_morse_potential_gradient(atoms, morse)
                istart, jstart = 3*i, 3*j
                istop, jstop = istart+3, jstart+3
                forces[istart:istop] -= g[0:3]
                forces[jstart:jstop] -= g[3:6]
            for bond in self.bonds:
                i, j, g = ff.get_bond_potential_gradient(atoms, bond)
                istart, jstart = 3*i, 3*j
                istop, jstop = istart+3, jstart+3
                forces[istart:istop] -= g[0:3]
                forces[jstart:jstop] -= g[3:6]
            for angle in self.angles:
                i, j, k, g = ff.get_angle_potential_gradient(atoms, angle)
                istart, jstart, kstart = 3*i, 3*j, 3*k
                istop, jstop, kstop = istart+3, jstart+3, kstart+3
                forces[istart:istop] -= g[0:3]
                forces[jstart:jstop] -= g[3:6]
                forces[kstart:kstop] -= g[6:9]
            for dihedral in self.dihedrals:
                i, j, k, l, g = ff.get_dihedral_potential_gradient(atoms, dihedral)
                istart, jstart, kstart, lstart = 3*i, 3*j, 3*k, 3*l
                istop, jstop, kstop, lstop = istart+3, jstart+3, kstart+3, lstart+3
                forces[istart:istop] -= g[0:3]
                forces[jstart:jstop] -= g[3:6]
                forces[kstart:kstop] -= g[6:9]
                forces[lstart:lstop] -= g[9:12]
            for vdw in self.vdws:
                i, j, g = ff.get_vdw_potential_gradient(atoms, vdw)
                istart, jstart = 3*i, 3*j
                istop, jstop = istart+3, jstart+3
                forces[istart:istop] -= g[0:3]
                forces[jstart:jstop] -= g[3:6]
            for coulomb in self.coulombs:
                i, j, g = ff.get_coulomb_potential_gradient(atoms, coulomb)
                istart, jstart = 3*i, 3*j
                istop, jstop = istart+3, jstart+3
                forces[istart:istop] -= g[0:3]
                forces[jstart:jstop] -= g[3:6]
            self.results['forces'] = np.reshape(forces, (len(atoms),3))
        if 'hessian' not in self.results:
            hessian = np.zeros((3*len(atoms),3*len(atoms)))
            for morse in self.morses:
                i, j, h = ff.get_morse_potential_hessian(atoms, morse)
                istart, jstart = 3*i, 3*j
                istop, jstop = istart+3, jstart+3
                hessian[istart:istop,istart:istop] += h[0:3,0:3]
                hessian[istart:istop,jstart:jstop] += h[0:3,3:6]
                hessian[jstart:jstop,istart:istop] += h[3:6,0:3]
                hessian[jstart:jstop,jstart:jstop] += h[3:6,3:6]
            for bond in self.bonds:
                i, j, h = ff.get_bond_potential_hessian(atoms, bond)
                istart, jstart = 3*i, 3*j
                istop, jstop = istart+3, jstart+3
                hessian[istart:istop,istart:istop] += h[0:3,0:3]
                hessian[istart:istop,jstart:jstop] += h[0:3,3:6]
                hessian[jstart:jstop,istart:istop] += h[3:6,0:3]
                hessian[jstart:jstop,jstart:jstop] += h[3:6,3:6]
            for angle in self.angles:
                i, j, k, h = ff.get_angle_potential_hessian(atoms, angle)
                istart, jstart, kstart = 3*i, 3*j, 3*k
                istop, jstop, kstop = istart+3, jstart+3, kstart+3
                hessian[istart:istop,istart:istop] += h[0:3,0:3]
                hessian[istart:istop,jstart:jstop] += h[0:3,3:6]
                hessian[istart:istop,kstart:kstop] += h[0:3,6:9]
                hessian[jstart:jstop,istart:istop] += h[3:6,0:3]
                hessian[jstart:jstop,jstart:jstop] += h[3:6,3:6]
                hessian[jstart:jstop,kstart:kstop] += h[3:6,6:9]
                hessian[kstart:kstop,istart:istop] += h[6:9,0:3]
                hessian[kstart:kstop,jstart:jstop] += h[6:9,3:6]
                hessian[kstart:kstop,kstart:kstop] += h[6:9,6:9]
            for dihedral in self.dihedrals:
                i, j, k, l, h = ff.get_dihedral_potential_hessian(atoms, dihedral)
                istart, jstart, kstart, lstart = 3*i, 3*j, 3*k, 3*l
                istop, jstop, kstop, lstop = istart+3, jstart+3, kstart+3, lstart+3
                hessian[istart:istop,istart:istop] += h[0:3,0:3]
                hessian[istart:istop,jstart:jstop] += h[0:3,3:6]
                hessian[istart:istop,kstart:kstop] += h[0:3,6:9]
                hessian[istart:istop,lstart:lstop] += h[0:3,9:12]
                hessian[jstart:jstop,istart:istop] += h[3:6,0:3]
                hessian[jstart:jstop,jstart:jstop] += h[3:6,3:6]
                hessian[jstart:jstop,kstart:kstop] += h[3:6,6:9]
                hessian[jstart:jstop,lstart:lstop] += h[3:6,9:12]
                hessian[kstart:kstop,istart:istop] += h[6:9,0:3]
                hessian[kstart:kstop,jstart:jstop] += h[6:9,3:6]
                hessian[kstart:kstop,kstart:kstop] += h[6:9,6:9]
                hessian[kstart:kstop,lstart:lstop] += h[6:9,9:12]
                hessian[lstart:lstop,istart:istop] += h[9:12,0:3]
                hessian[lstart:lstop,jstart:jstop] += h[9:12,3:6]
                hessian[lstart:lstop,kstart:kstop] += h[9:12,6:9]
                hessian[lstart:lstop,lstart:lstop] += h[9:12,9:12]
            for vdw in self.vdws:
                i, j, h = ff.get_vdw_potential_hessian(atoms, vdw)
                istart, jstart = 3*i, 3*j
                istop, jstop = istart+3, jstart+3
                hessian[istart:istop,istart:istop] += h[0:3,0:3]
                hessian[istart:istop,jstart:jstop] += h[0:3,3:6]
                hessian[jstart:jstop,istart:istop] += h[3:6,0:3]
                hessian[jstart:jstop,jstart:jstop] += h[3:6,3:6]
            for coulomb in self.coulombs:
                i, j, h = ff.get_coulomb_potential_hessian(atoms, coulomb)
                istart, jstart = 3*i, 3*j
                istop, jstop = istart+3, jstart+3
                hessian[istart:istop,istart:istop] += h[0:3,0:3]
                hessian[istart:istop,jstart:jstop] += h[0:3,3:6]
                hessian[jstart:jstop,istart:istop] += h[3:6,0:3]
                hessian[jstart:jstop,jstart:jstop] += h[3:6,3:6]
            self.results['hessian'] = hessian 

    def get_hessian(self, atoms=None):
        return self.get_property('hessian', atoms)
