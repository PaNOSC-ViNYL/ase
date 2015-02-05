.. highlight:: bash

.. index:: command line tools

.. _cli:

=========================
Command line tools in ASE
=========================

ASE has the following command line tools:
    
* :ref:`ase-gui`: graphical user interface
* :ref:`ase-db`: manipulation of databases
* ase-build: build simple molecule or bulk structure
* ase-run: run calculations with ASE's calculators
* ase-info


Help
====

For all command-line tools, you can do::
    
    ase-xxx --help
    
or ``ase-xxx -h`` to get help.


.. _bash completion:
    
Bash completion
===============

You can enable bash completion by adding this line to your ``~/.bashrc``::
    
    complete -o default -C _ase_bash_complete.py ase-db ase-run ase-build ase-info ase-gui
