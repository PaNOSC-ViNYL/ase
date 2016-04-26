.. highlight:: bash

.. index:: command line tools

.. _cli:

==================
Command line tools
==================

ASE has the following command line tools:
    
* :ref:`ase-gui`: graphical user interface
* :ref:`ase-db`: manipulation of databases
* ase-build: build simple molecule or bulk structure
* ase-run: run calculations with ASE's calculators
* ase-info


Python -m tricks
================

Some ASE modules can be invoked directly form the command line using ``python
-m``.
    
:ref:`stylecheck`::
    
    $ python -m ase.utils.stylecheck source.py

Equation of state::
    
    $ python -m ase.eos [-p] traj-file, ...
    
:ref:`iso surface`::

    $ python -m ase.visulaize.mlab [options] filename
    
Determine file type(s)::
    
    $ python -m ase.io.formats file ...

Convert old db-files to new::
    
    $ python -m ase.db.convert db-file
    
Show content of aff-file::
    
    $ python -m ase.io.aff [options] aff-file [item number]
    
:ref:`convert`::
    
    $ python -m ase.io.pickletrajectory a1.traj [a2.traj ...]


Help
====

For all command-line tools, you can do::
    
    $ ase-gui --help
    $ python -m ase.eos --help
    
to get help (or ``-h`` for short).


.. _bash completion:
    
Bash completion
===============

You can enable bash completion by adding this line to your ``~/.bashrc``::
    
    complete -o default -C _ase_bash_complete.py ase-db ase-run ase-build ase-info ase-gui
