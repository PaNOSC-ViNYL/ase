from collections import defaultdict

import matplotlib.pyplot as plt

import ase.db


energies = defaultdict(list)
for row in ase.db.connect('results.db').select():
    energies[row.formula].append(row.energy)
emin = {formula: min(energies[formula]) for formula in energies}

n = defaultdict(list)
for row in ase.db.connect('results.db').select(sort='sid'):
    if row.energy - emin[row.formula] < 0.01:
        nsteps = row.n
    else:
        nsteps = float('inf')
    n[row.optimizer].append(nsteps)

N = sorted(n.items(), key=lambda x: sum(x[1]))
for o, n in N:
    print('{:18}{}'.format(o, ' '.join('{:3}'.format(x) for x in n)))

plt.imshow([n for o, n in N])
plt.show()
