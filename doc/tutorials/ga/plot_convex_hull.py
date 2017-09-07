from ase.phasediagram import PhaseDiagram
from ase.db import connect

db = connect('hull.db')

refs = []
for dct in db.select('relaxed=1'):
    refs.append((dct.formula, dct.raw_score))

pd = PhaseDiagram(refs)
pd.plot()
