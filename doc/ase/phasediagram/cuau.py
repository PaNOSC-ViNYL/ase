# creates: cuau.svg
from ase.phasediagram import PhaseDiagram
refs = [('Cu', 0.0),
        ('Au', 0.0),
        ('CuAu', -0.5),
        ('Cu2Au', -0.7),
        ('CuAu2', -0.2)]
pd = PhaseDiagram(refs)
pd.plot()
import matplotlib.pyplot as plt
plt.savefig('cuau.svg')
print(pd.decompose('Cu3Au'))
