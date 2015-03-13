from __future__ import print_function
from ase.test import NotAvailable

raise NotAvailable

import numpy as np

import ase.db
from ase.phasediagram import bisect, Pourbaix

if 0:
    N = 80
    A = np.zeros((N, N), int)
    A[:] = -1
    
    def f(x, y):
        dmin = 100
        for i, (a, b) in enumerate([(0, 0), (0, 2), (1, 1)]):
            d = (x - a)**2 + (y - b)**2
            if d < dmin:
                dmin = d
                imin = i
        return imin

    bisect(A, np.linspace(0, 2, N), np.linspace(0, 2, N), f)
    print(A)

    import matplotlib.pyplot as plt
    plt.imshow(A)
    plt.show()
    

if 0:
    con = ase.db.connect('cubic_perovskites.db')
    references = [(row.count_atoms(), row.energy)
                  for row in con.select('reference')]
    std = {}
    for count, energy in references:
        if len(count) == 1:
            symbol, n = list(count.items())[0]
            assert symbol not in std
            std[symbol] = energy / n

    std['O'] += 2.46

    refs = []
    for refcount, energy in references:
        for symbol, n in refcount.items():
            energy -= n * std[symbol]
        if list(refcount) == ['O']:
            energy = 0.0
        refs.append((refcount, energy))
else:
    refs = [({'O': 1}, 0.0),
            ({'O': 4, 'Ti': 2}, -17.511826939900217),
            ({'Sr': 4, 'O': 4}, -20.474907588620653),
            ({'Sr': 4}, 0.0),
            ({'Ti': 2}, 0.0)]

pb = Pourbaix(refs, Sr=1, Ti=1, O=3)
print(pb.decompose(0, 0))
pH = np.linspace(-1, 15, 17)
if 0:
    d, names = pb.diagram([0], pH)
    print(d)
    print('\n'.join(names))
U = np.linspace(-2, 2, 5)
if 0:
    d, names = pb.diagram(U, [0])
    for i, u in zip(d, U):
        print(u, names[i])

if 0:
    U = np.linspace(-3, 3, 200)
    pH = np.linspace(-1, 15, 300)
    d, names = pb.diagram(U, pH, plot=True)
