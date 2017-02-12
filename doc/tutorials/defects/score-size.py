# creates: score-size-fcc2sc.svg

import matplotlib.pyplot as plt
import json
import os
import glob

for fname in glob.glob('Popt-*.json'):
    tag = os.path.basename(fname).replace('Popt-', '').replace('.json', '')

    with open(fname) as data_file:
        data = json.load(data_file)
    x = []
    y = []
    for nuc, rec in sorted(data.items()):
        x.append(nuc)
        y.append(rec['dev'])

    plt.figure(figsize=(4, 3))
    plt.text(1950, 0.6,
             tag.replace('2', r' $\rightarrow$ '), horizontalalignment='right')
    plt.xlabel(r'Number of primitive unit cells $N_{uc}$')
    plt.ylabel(r'Optimality measure $\bar \Delta$')
    plt.axis([0, 2000, -0.05, 0.7])
    plt.plot(x, y, 'bo')
    plt.savefig('score-size-%s.svg' % tag, bbox_inches='tight')
