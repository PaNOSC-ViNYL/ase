The Atoms object
================

.. automodule:: ase.atoms

.. autoclass:: Atoms

Working with the methods of :class:`~ase.atoms.Atoms`
---------------------------------------------------------

Like with a single :class:`~ase.atoms.Atom` the properties of a collection of atoms
can be accessed and changed with get- and set-methods. For example
the positions of the atoms can be addressed as

>>> a = Atoms('N3', [(0, 0, 0), (1, 0, 0), (0, 0, 1)])
>>> a.get_positions()
array([[ 0.,  0.,  0.],
       [ 1.,  0.,  0.],
       [ 0.,  0.,  1.]])
>>> a.set_positions([(2, 0, 0), (0, 2, 2), (2, 2, 0)])
>>> a.get_positions()
array([[ 2.,  0.,  0.],
       [ 0.,  2.,  2.],
       [ 2.,  2.,  0.]])

It is also possible to work directly with the attribute "position" for
a single atom or "positions" for several atoms. Here we change the
position of the 2nd atom (which has count number 1 because Python
starts at zero):

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
default is the 3x3 unit matrix as can be seen from

>>> a.get_cell()
array([[ 1.,  0.,  0.],
       [ 0.,  1.,  0.],
       [ 0.,  0.,  1.]])


The cell can be defined or changed using the
:epydoc:`ase.atoms.Atoms.set_cell` method. Changing the unit cell
does per default not move the atoms:

>>> a.set_cell(2 * identity(3))
>>> a.get_cell()
array([[ 2.,  0.,  0.],
       [ 0.,  2.,  0.],
       [ 0.,  0.,  2.]])
>>> a.get_positions()
array([[ 2.,  0.,  0.],
       [ 1.,  1.,  0.],
       [ 2.,  2.,  0.]])

However if we set fix=False the atomic positions are scaled with the unit cell:

>>> a.set_cell(identity(3), fix=False)
>>> a.get_positions()
array([[ 1. ,  0. ,  0. ],
       [ 0.5,  0.5,  0. ],
       [ 1. ,  1. ,  0. ]])

The :epydoc:`ase.atoms.Atoms.set_pbc` method specifies whether
periodic boundary conditions are to be used in the directions of the
three vectors of the unit cell. With the default unit cell
[(1,0,0),(0,1,0),(0,0,1)] a slab calculation with periodic boundary
conditions in x and y and free boundary condtions in z is obatined
through

>>> a.set_pbc((True, True, False))

A calculator can be attached to the atoms with the purpose
of calculating energies and forces on the atoms. ASE works with many
different :mod:`calculators`.

A calculator object *calc* is attached to the atoms like this:

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

In case of a DFT calculator, it is up to the user to check exactly what
the :epydoc:`ase.atoms.Atoms.get_potential_energy` method returns. For
example it may be the result of a calculation with a finite
temperature smearing of the occupation numbers extrapolated to zero
temperature. More about this can be found for the different
:mod:`calculators` XXX Is get_potential_energy well defined for the
different calculators ? XXX

More information about how to manipulate 
