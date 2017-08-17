.. highlight:: bash

.. index:: command line tools

.. _cli:

=================
Command line tool
=================

ASE has a command line tool called :program:`ase` with the following
sub-commands:

==============  =================================================
sub-command     description
==============  =================================================
help            Help for sub-command
info            Print information about files or system
test            Test ASE
gui             ASE's :ref:`graphical user interface <ase-gui>`
find            Find files with atoms in
db              Manipulate and query :ref:`ASE database <ase-db>`
run             Run calculation with one of ASE's calculators
build           Build an atom, molecule or bulk structure
eos             Calculate equation of state
ulm             Show content of ulm-file
nomad-upload    Upload files to NOMAD
band-structure  Plot band-structure
completion      Add tab-completion for Bash
==============  =================================================


Python -m tricks
================

Some ASE modules can be invoked directly form the command line using ``python3
-m``.

:ref:`stylecheck`::

    $ python -m ase.utils.stylecheck source.py

:ref:`iso surface`::

    $ python -m ase.visulaize.mlab [options] filename

Convert old db-files to new::

    $ python -m ase.db.convert db-file

:ref:`convert`::

    $ python -m ase.io.pickletrajectory a1.traj [a2.traj ...]


Help
====

For all command-line tools, you can do::

    $ ase --help
    $ ase sub-command --help
    $ python -m module --help

to get help (or ``-h`` for short).


.. _bash completion:

Bash completion
===============

You can enable bash completion like this::

    $ ase completions

This will append a line like this::

    complete -o default -C /path/to/ase/ase/cli/complete.py ase

to your ``~/.bashrc``.
