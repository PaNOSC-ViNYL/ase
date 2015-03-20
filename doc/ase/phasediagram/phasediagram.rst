.. module:: ase.phasediagram

====================================
Phase diagrams and Pourbaix diagrams
====================================

.. autoclass:: ase.phasediagram.PhaseDiagram

Here is a simple example using some made up numbers for Cu-Au alloys:
    
>>> refs = [('Cu', 0.0),
...         ('Au', 0.0),
...         ('CuAu', -0.5),
...         ('Cu2Au', -0.7),
...         ('Cu2Au', -0.2)]
>>> pd = PhaseDiagram(refs)
Species: Au, Cu
References: 5
0    Cu             0.000
1    Au             0.000
2    CuAu          -0.500
3    Cu2Au         -0.700
4    CuAu2         -0.200
Simplices: 3

The phase diagram looks like this:
    
>>> pd.plot()

.. image:: cuau.png

.. automethod:: PhaseDiagram.plot

If you want to see what CuAu will decompose into, you can use the
:meth:`~PhaseDiagram.decompose` method:
    
.. automethod:: PhaseDiagram.decompose

>>> pd.decompose('Cu3Au')  # or pd.decompose(Cu=3, Au=1)
reference    coefficient      energy
------------------------------------
Cu                     1       0.000
Cu2Au                  1      -0.700
------------------------------------
Total energy:                 -0.700
------------------------------------
(-0.69999999999999996, array([0, 3], dtype=int32), array([ 1.,  1.]))


Pourbaix diagrams
=================

Let's create a Pourbaix diagram for ZnO from experimental numbers.

>>> from ase.phasediagram import Pourbaix, solvated
>>> refs = solvated('Zn')
>>>> print(refs)
[('HO2Zn-(aq)', -4.801274772854441), ('O2Zn--(aq)', -4.0454382546928365), ('HOZn+(aq)', -3.5207324675582736), ('OZn(aq)', -2.9236086089762137), ('H2O(aq)', -2.458311658897383), ('Zn++(aq)', -1.5264168353005447), ('H+(aq)', 0.0)]

We use the :func:`solvated` function to get solvation energies for zinc
containing molecules (plus water and a proton):
    
.. autofunction:: solvated

We add two solids and one more disolved molecule to the references and create
a :class:`Pourbaix` object:
    
>>> refs += [('Zn', 0.0), ('ZnO', -3.323), ('ZnO2(aq)', -2.921)]
>>> pb = Pourbaix(refs, Zn=1, O=1)

To see what ZnO will :meth:`~Pourbaix.decompose` to at a potential of 1 eV
and a pH of 9.0, we do this:
    
>>> coefs, energy = pb.decompose(1.0, 9.0)
0    HO2Zn-(aq)    -5.158
1    O2Zn--(aq)    -4.403
2    HOZn+(aq)     -3.878
3    OZn(aq)       -3.281
4    H2O(aq)       -2.458
5    Zn++(aq)      -1.884
6    H+(aq)        -0.536
7    Zn             0.000
8    OZn           -3.323
9    O2Zn(aq)      -3.278
10   e-            -1.000
reference    coefficient      energy
------------------------------------
H2O(aq)               -1      -2.458
H+(aq)                 2      -0.536
O2Zn(aq)               1      -3.278
e-                     2      -1.000
------------------------------------
Total energy:                 -3.891
------------------------------------
>>> print(coefs, energy)
(array([  0.00000000e+00,   0.00000000e+00,   6.66133815e-16,
         0.00000000e+00,  -1.00000000e+00,   0.00000000e+00,
         2.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         1.00000000e+00,   2.00000000e+00]), -3.8913313372636829)

The full :meth:`~Pourbaix.diagram` is calculated like this:
    
>>> import numpy as np
>>> U = np.linspace(-2, 2, 200)
>>> pH = np.linspace(-2, 16, 300)
>>> d, names, text = pb.diagram(U, pH, plot=True)

.. image:: zno.png

.. autoclass:: ase.phasediagram.Pourbaix
    :members:
    :member-order: bysource
