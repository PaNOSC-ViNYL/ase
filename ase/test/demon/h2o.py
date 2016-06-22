import ase.calculators.demon as demon
from ase import Atoms
from ase.optimize import BFGS
import numpy as np

tol = 1.0e-8

# d = 0.9575
d = 0.9775
# t = np.pi / 180 * 104.51
t = np.pi / 180 * 110.51
atoms = Atoms('H2O',
              positions=[(d, 0, 0),
                         (d * np.cos(t), d * np.sin(t), 0),
                         (0, 0, 0)])

# set up deMon calculator
basis = {'all': 'aug-cc-pvdz',
         'O': 'RECP6|SD'}
auxis = {'all': 'GEN-A2*'}
input_arguments = {'GRID': 'FINE'}
    
calc = demon.Demon(basis=basis,
                   auxis=auxis,
                   scftype='RKS TOL=1.0E-6 CDF=1.0E-5',
                   guess='TB',
                   xc=['BLYP', 'BASIS'],
                   input_arguments=input_arguments)

atoms.set_calculator(calc)

# energy
energy = atoms.get_potential_energy()

ref = -469.604903217
print('energy')
print(energy)
error = np.sqrt(np.sum((energy - ref)**2))
print('diff from reference:')
print(error)

assert(error < tol)

# dipole
dipole = atoms.get_dipole_moment()

ref = np.array([0.19228183, 0.27726241, 0.0])
error = np.sqrt(np.sum((dipole - ref)**2))
print('dipole')
print(dipole)
print('diff from reference:')
print(error)

assert(error < tol)


# numerical forces
forces_num = calc.calculate_numerical_forces(atoms, d=0.001)

ref = np.array([[-1.26056790e-01, 4.10007704e-01, 2.85719636e-04],
                [4.28062465e-01, 2.56059233e-02, 2.17691195e-04],
                [-3.02019280e-01, -4.35613627e-01, -5.03410803e-04]])
error = np.sqrt(np.sum((forces_num - ref)**2))
print('forces_num')
print(forces_num)
print('diff from reference:')
print(error)

assert(error < tol)


# analytical forces
forces_an = atoms.get_forces()

ref = np.array([[-1.26446897e-01, 4.09628295e-01, -0.00000000e+00],
                [4.27934556e-01, 2.50425533e-02, -5.14220807e-05],
                [-2.99225088e-01, -4.31534101e-01, -5.14220807e-05]])

error = np.sqrt(np.sum((forces_an - ref)**2))
print('forces_an')
print(forces_an)
print('diff from reference:')
print(error)

assert(error < tol)

# optimize geometry
dyn = BFGS(atoms)
dyn.run(fmax=0.01)

positions = atoms.get_positions()
ref = np.array([[9.61364575e-01, 2.81689442e-02, -1.58730812e-06],
                [-3.10444390e-01, 9.10289259e-01, -5.66399225e-06],
                [-1.56957804e-02, -2.26044113e-02, -2.34155677e-06]])
error = np.sqrt(np.sum((positions - ref)**2))
print('positions')
print(positions)
print('diff from reference:')
print(error)

assert(error < tol)

print('tests passed')



