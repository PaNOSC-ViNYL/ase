# creates: score-size-fcc.svg

import matplotlib.pyplot as plt
import json
import glob
import os

for fname in glob.glob('Popt*.json'):
    tag = os.path.basename(fname).replace('Popt-', '').replace('.json','')

    with open(fname) as data_file:    
        data = json.load(data_file)
    x = []
    y = []
    for nuc,rec in sorted(data.items()):
        x.append(nuc)
        y.append(rec['score'])

    plt.figure(figsize=(5, 5))
    plt.xlabel('Number of primitive unit cells')
    plt.ylabel('Score $\Delta$')
    plt.plot(x, y, 'bo')
    plt.savefig('score-size-%s.pdf' % tag, bbox_inches='tight')
