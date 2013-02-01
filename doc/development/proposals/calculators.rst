=============================
Calculator interface proposal
=============================

All ASE calculators should behave similarly if there is no good reason
for them not to.  This should make it simpler for both users and developers.

This proposal tries to define how a good ASE calculator should behave.
The goal is to have ASE calculators:

* that share more code
* are more uniform to use
* are better tested

Setting some standards is a good thing, but we should also be careful
not to set too strict rules that could limit each calculator to the
lowest common denominator.

So far, this proposal is mostly ideas and questions.  


Behavior
========

When a calculator calculates the energy, forces, stress tensor, total
magnetic momen, atomic magnetic moments or dipole moment, it should
store a copy of the system (atomic numbers, atomic positions, unit
cell and boundary conditions).  When asked again, it should return the
value already calculated if the system hasn't been changed.

If calculational parameters such as plane wave cutoff or XC-functional
has been changed the calculator should throw away old calculated
values.


Standards parameters
====================

The standard keywords that all calculators must use (if they make
sense) are: ``xc``, ``kpts``, ``smearing``, ``width``, ``charge`` and
``nbands``.  Each calculator will have its own default values for
these parameters --- see recommendations below.  I addition,
calculators will typically have many other parameters.  The units are
eV and Å.

Initial magnetic moments are taken from the :class:`~ase.atoms.Atoms`
object.

:xc:

  It is recommended that ``'LDA'`` and ``'PBE'`` are valid options.

:kpts:

  * ``None``: Gamma-point
  
  * ``(n1,n2,n3)``: Monkhorst-Pack grid
  
  * ``(n1,n2,n3,'gamma')``: Shifted Monkhorst-Pack grid that includes `\Gamma`
  
  * ``[(k11,k12,k13),(k21,k22,k23),...]``: Explicit list in units of the
    reciprocal lattice vectors
  
  * ``kpts=3.5``: `\vec k`-point density as in 3.5 `\vec k`-points per
    Å\ `^{-1}`.

:smearing:

  The smearing parameter can be one of these string:

  * ``'Fermi-Dirac'`` or ``'FD'``
  * ``'Gaussian'``
  * ``'Methfessel-Paxton-n'`` or ``'MPn'``, where `n` is the order
    (`n=0` is the same as ``'Gaussian'``)
  * or lower-case versions of any of the above

  Recommended value: ``'fermi-dirac'`` with ``width=0.1`` eV.

:width:

  The width parameter used for the chosen smearing method (in eV).
  Recommended value is 0.1 eV combined with ``smearing='fermi-dirac'``

:charge:

  Charge of the system in units of `|e|` (``charge=1`` means one
  electron has been removed).


:nbands:

  Each band can be occupied by two electrons.

  
ABC calculator example
======================

A calculator should be able to prefix all output files with a given
label or run the calculation in a directory with a specified name.
There are three possibilities for the first argument (called
``label``) of the constructor of a calculator object:

* Name of a file containing all results of a calculation (possibly
  containing a directory).

* A prefix used for several files containing results.  The label may
  have both a directory part and a prefix part like ``'LDA/mol1'``.

* Name of a directory containing result files with fixed names.

Each calculator can decide what the default value is: ``None`` for no
output, ``'-'`` for standard output or something else.  All others
parameters are given as keyword arguments.

Example:  Do a calculation with ABC calculator and write results to
:file:`si.abc`:

>>> atoms = ...
>>> atoms.calc = ABC('si.abc', xc='LDA')
>>> atoms.get_potential_energy()
-1.2

The default behavior of of reading from :file:`si.abc` and also
writing results from following calculations to the same file can be
changed by using the ``output`` keyword (``output='si-new.abc'``).

An alternative way to connect atoms and calculator:

>>> atoms = ...
>>> calc = ABC('si.abc', atoms=atoms, xc='LDA')
>>> atoms.get_potential_energy()
-1.2

This will automatically attach the calculator to the atoms and if the
file :file:`si.abc` exists the atoms will be updated form that file.
This will allow you to use the same script to do the initial
calculation where :file:`si.abc` does not exist and following
calculations where atoms may have been moved arround be an
optimization algorithm.

Read atoms with ABC calculator attaced from a previous calculation:

>>> atoms = ABC.read_atoms('si.abc')
>>> atoms.calc
<ABC-calculator>
>>> atoms.get_potential_energy()
-1.2

The class method :meth:`read_atoms()` is equivalent to:

>>> atoms = ABC('si.abc').get_atoms()

If we do:

>>> atoms = ABC.read('si.abc')
>>> atoms.rattle()            # change positions and/or
>>> atoms.calc.set(xc='PBE')  # change a calculator-parameter
>>> atoms.get_potential_energy()
-0.7

then the :file:`si.abc` will be overwritten or maybe appended to.

The command used to start the ABC code must be given in an environment
variable called :envvar:`ASE_ABC_COMMAND` or as a ``command``
keyword.  The command can look like this::

  mpiexec abc LABEL.input > LABEL.output

or like this::

  ~/bin/start_abc.py LABEL

The ``LABEL`` strings will be substituted by the ``label`` prefix.
If neither the environment variable or the ``command`` keyword is
specified, the calculator will raise a ``NotAvailable`` exception,
which will make the test-suite skip such tests.


Pre- and post-run hooks
=======================

Calculators can call call-back functions before and after a
calculation.  These are stored in a dictionary called ``callbacks``::

  {'before': [], 'after': []}

The two lists can contain any number of functions specified as
``(function, args, kwargs)`` tuples.  The lists can be manipulated
directly or apended to by using the ``attach_callback(when, function,
*args, **kwargs)`` method, where ``when`` is one of ``'before'`` or
``'after'``.


Implementation
==============

* Common base class for all calculators: ``Calculator``.  Takes care
  of restart from file logic, handles setting of parameters and checks
  for state changes.

* A ``FileIOCalculator`` for the case where we need to:

  * write input file(s)
  * run Fortran code
  * read output file(s)

* Helper function to deal with ``kpts`` keyword.
