import numpy as np
from ase.phasediagram import PhaseDiagram
from ase.db import connect
from ase.visualize import view

db = connect('hull.db')

refs = []
dcts = list(db.select('relaxed=1'))
for dct in dcts:
    refs.append((dct.formula, dct.raw_score))

pd = PhaseDiagram(refs)
pd.plot(only_label_simplices=True)

# View the simplices of the convex hull
simplices = []
toview = sorted(np.array(dcts)[pd.hull], key=lambda x: x.mass)
for dct in toview:
    simplices.append(dct.toatoms())

view(simplices)
