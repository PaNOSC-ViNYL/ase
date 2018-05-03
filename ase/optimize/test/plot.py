from collections import defaultdict
from math import inf

import ase.db


def analyze(filename, tag='results'):
    energies = defaultdict(list)
    formulas = []
    db = ase.db.connect(filename)
    for row in db.select(sort='name'):
        if row.formula not in formulas:
            formulas.append(row.formula)
        energies[row.formula].append(row.get('energy', inf))
    emin = {formula: min(energies[formula]) for formula in energies}

    data = defaultdict(list)
    for row in db.select(sort='name'):
        if row.get('energy', inf) - emin[row.formula] < 0.01:
            nsteps = row.n if row.n < 100 else 9999
            t = row.t
        else:
            nsteps = 9999
            t = 9999
        data[row.optimizer].append((nsteps, t))

    print(formulas)

    D = sorted(data.items(), key=lambda x: sum(y[0] for y in x[1]))
    with open(tag + '-iterations.csv', 'w') as f:
        print('optimizer,' + ','.join(formulas), file=f)
        for o, d in D:
            print('{:18},{}'
                  .format(o, ','.join('{:3}'.format(x[0])
                                      if x[0] < 100 else '   '
                                      for x in d)),
                  file=f)

    D = sorted(data.items(), key=lambda x: sum(y[1] for y in x[1]))
    with open(tag + '-time.csv', 'w') as f:
        print('optimizer,' + ','.join(formulas), file=f)
        for o, d in D:
            print('{:18},{}'
                  .format(o, ','.join('{:8.3f}'.format(x[1])
                                      if x[0] < 100 else '        '
                                      for x in d)),
                  file=f)


if __name__ == '__main__':
    analyze('results.db')
