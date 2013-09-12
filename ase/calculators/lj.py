import numpy as np
from ase.calculators.calculator import Calculator


class LennardJones(Calculator):
    implemented_properties = ['energy', 'forces']
    default_parameters = {'epsilon': 1.0,
                          'sigma': 1.0}
    nolabel = True

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

    def calculate(self, atoms, properties, changes):
        epsilon = self.parameters.epsilon
        sigma = self.parameters.sigma
        sigma2 = sigma**2
        
        positions = atoms.get_positions()
        energy = 0.0
        forces = np.zeros((len(atoms), 3))

        for i in range(len(positions) - 1):
            dist = positions[i, :] - positions[i + 1:, :]
            d2 = (dist**2).sum(axis=1)
            c6 = (sigma2 / d2) ** 3
            c12 = c6 ** 2

            energy += sum(4 * epsilon * (c12 - c6))
            ljf = (24 * epsilon * (2 * c12 - c6) / d2)
            ljf = ljf.reshape((len(d2), 1)) * dist
            
            forces[i + 1:, :] -= ljf
            forces[i, :] += ljf.sum(axis=0)

        self.results['energy'] = energy
        self.results['forces'] = forces
