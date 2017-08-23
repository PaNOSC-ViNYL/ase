""" Water dimer calculation in which each molecule is calculated quantum 
mechanically and the interaction between the molecules is electrostatic. The process is repeated until self consitence. """

from numpy.linalg import norm
from ase.collections import s22
from ase.calculators.turbomole import Turbomole

params = {'esp fit': 'kollman', 'multiplicity': 1}
dimer = s22['Water_dimer']

# system partitioning
part1 = dimer[0:3]
part2 = dimer[3:6]

def polarization_cycle(part1, part2, c2=None):
    prop = {}
    calc = Turbomole(atoms=part1, **params)
    if c2 is not None:
        calc.embed(charges=c2, positions=part2.positions)
    prop['e1'] = part1.get_potential_energy()
    prop['c1'] = part1.get_charges()
    calc = Turbomole(atoms=part2, **params)
    calc.embed(charges=prop['c1'], positions=part1.positions)
    prop['e2'] = part2.get_potential_energy()
    prop['c2'] = part2.get_charges()
    return prop

new = polarization_cycle(part1, part2)
prop = {'e1': [], 'e2': [], 'c1': [], 'c2': []}
for key in prop:
    prop[key].append(new[key])

# start values and convergence criteria
conv = {'e': [1.0], 'c': [1.0]}
thres = {'e': 1e-4, 'c': 1e-2}
iteration = 0
while any([conv[key][-1] > thres[key] for key in conv]):
    iteration += 1
    new = polarization_cycle(part1, part2, c2=prop['c2'][-1])
    for key in prop:
        prop[key].append(new[key])

    conv['e'].append(
        (abs(prop['e1'][-1]-prop['e1'][-2])
        +abs(prop['e2'][-1]-prop['e2'][-2]))/
        (abs(prop['e1'][-2])
        +abs(prop['e2'][-2]))
    )
    conv['c'].append(
        (norm(prop['c1'][-1]-prop['c1'][-2])
        +norm(prop['c2'][-1]-prop['c2'][-2]))/
        (norm(prop['c1'][-2])
        +norm(prop['c2'][-2]))
    )
    fmt = 'iteration {0:d}: convergence of energy {1:10e}; of charge {2:10e}'
    print(fmt.format(iteration, conv['e'][-1], norm(conv['c'][-1])))

# check the result
ref = {
    'e1': -2077.7082947500003,
    'e2': -2077.3347674372353,
    'c1': [-0.133033,  0.238218, -0.105186],
    'c2': [-0.844336,  0.422151,  0.422184]
}

dev = {}
for key in ref:
    val = prop[key][-1]
    measurement = val if isinstance(val, float) else norm(val)
    reference = ref[key] if isinstance(ref[key], float) else norm(ref[key])
    dev[key] = (measurement-reference)/reference
    print('Deviation of {0} is {1:10f}'.format(key, dev[key]))

# allow deviations of up to 5%
assert all([dev[key] < 5e-2 for key in dev]), 'deviation too large'
