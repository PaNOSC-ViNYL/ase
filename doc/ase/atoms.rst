The Atoms object
================

.. automodule:: ase.atoms

.. autoclass:: Atoms

Working with the methods of :class:`~ase.atoms.Atoms`
---------------------------------------------------------

Like with a single :class:`~ase.atoms.Atom` the properties of a collection of atoms
can be accessed and changed with "get-" and "set-"methods. For example
the positions of the atoms can be addressed as

>>> a=Atoms('N3',[(0,0,0),(1,0,0),(0,0,1)])
>>> a.get_positions()
array([[ 0.,  0.,  0.],
       [ 1.,  0.,  0.],
       [0., 0., 1].])
>>> a.set_positions([(2,0,0),(0,2,2),(2,2,0)])
>>> a.get_positions()
array([[ 2.,  0.,  0.],
       [ 0.,  2.,  2.],
       [2., 2., 0.]])

It is also possible to work directly with the attribute "position" for
a single atom or "positions" for several atoms. Here we change the
position of the 2nd atom (which has count number 1 because python
starts with zero):

>>> a[1].position=(1,1,0)
>>> a.get_positions()
array([[2., 0., 0.],
	    [1., 1., 0.],
	    [2., 2., 0.]])
>>> a.positions
array([[2., 0., 0.],
	    [1., 1., 0.],
	    [2., 2., 0.]])

The :class:`~ase.atoms.Atoms` object holds a unit cell which by
default is the 3x3 unit matrix. The cell can be set using the
:epydoc:`ase.atoms.Atoms.set_cell` method

The :epydoc:`ase.atoms.Atoms.set_pbc` method specifies whether periodic
boundary conditions are to be used in the directions of the three
vectors of the unit cell. With the default unit cell
[(1,0,0),(0,1,0),(0,0,1)] a slab calculation with periodic boundary
conditions in x and y and free boundary condtions in z is obatined
through

>>> a.set_pbc((True,True,False))

A calculator can be attached to the atoms with the purpose
of calculating energies and forces on the atoms. ASE works with many
different :mod:`calculators`.

A calculator object "calc" is attached to the atoms like this:

>>> a.set_calculator(calc)

After the calculator has been appropriately setup the energy of the
atoms can be obtained through

>>> a.get_potential_energy()

The term "potential energy" here means for example the total energy of
a DFT calculation, which includes both kinetic, electrostatic, and
exchange-correlation energy for the electrons. The reason it is called
potential energy is that the atoms might also have a kinetic energy
(from the moving nuclei) and that is obtained with

>>> a.get_kinetic_energy()

In case of a DFT calculator it is up to the user to check exactly what
the :epydoc:`ase.atoms.Atoms.get_potential_energy` method returns. For
example it may be the result of a calculation with a finite
temperature smearing of the occupation numbers extrapolated to zero
temperature. More about this can be found for the different
:mod:`calculators` XXX Is get_potential_energy well defined for the
different calculators ? XXX

More examples of manipulating atomic positions
--------------------------------------------------

We will end up with a one layer slab with one adatom

Define the slab atoms:

>>> from ase import *
>>> atoms = Atoms([Atom('Ni', (0, 0, 0)),
...                      Atom('Ni', (0.45, 0, 0)),
...                      Atom('Ni', (0, 0.5, 0)),
...                      Atom('Ni', (0.5, 0.5, 0))])


Have a look at the individual atoms:

>>> atoms[0]
Atom('Ni', [0.0, 0.0, 0.0], atoms=..., index=0)
>>> atoms[1]
Atom('Ni', [0.45, 0.0, 0.0], atoms=..., index=1)
>>> atoms[2]
Atom('Ni', [0.0, 0.5, 0.0], atoms=..., index=2)
>>> atoms[3]
Atom('Ni', [0.5, 0.5, 0.0], atoms=..., index=3)

Let us assume we forgot how many atoms we set up:

>>> atoms[4]
Traceback (most recent call last):
File "<stdin>", line 1, in ?
IndexError: list index out of range

Wrong because we only have four atoms

>>> len(atoms)
4

Change the position of the 2nd atom in the list

>>> atoms[1].set_position((0.5,0,0))
>>> atoms.get_positions()
array([[ 0. ,  0. ,  0. ],
       [ 0.5,  0. ,  0. ],
       [ 0. ,  0.5,  0. ],
       [ 0.5,  0.5,  0. ]])

What is the unit cell so far?

>>> atoms.get_cell()
array([[ 1.,  0.,  0.],
       [ 0.,  1.,  0.],
       [ 0.,  0.,  1.]])

Now, setup a p(2x2) cell in a hexagonal surface.
"a" is the fcc lattice constant, the cell is 10 layers high:





