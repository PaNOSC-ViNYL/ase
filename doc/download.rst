========================
Download and install ASE
========================

.. contents



Requirements for ASE
====================

The following packages are required for basic ASE functionality:

1) Python_.
2) NumPy_.

.. _Python: http://www.python.org
.. _NumPy: http://www.scipy.org/NumPy


It is highly recommended (but not required) to install also Matplotlib_
and pygtk_.  These two may already be installed on your system.


.. _Matplotlib: http://matplotlib.sourceforge.net
.. _pygtk: http://www.pygtk.org


Installation
============

.. highlight:: bash

Get the source code from svn::

  $ cd
  $ svn checkout https://svn.fysik.dtu.dk/projects/ase/trunk ase
  $ cd ase
	
or from this tarfile (python-ase-3.0.0.tar.gz_)::

  $ cd
  $ tar xtzf python-ase-3.0.0.tar.gz
  $ cd python-ase-3.0.0

.. _python-ase-3.0.0.tar.gz: python-ase-3.0.0.tar.gz

If you have root-access, you do this::

  $ python setup.py install

If you don't have root-access, you must put the directory
:file:`$HOME/ase3k` in your :envvar:`PYTHONPATH` environment variable,
and the directory :file:`$HOME/ase3k/tools` in your :envvar:`PATH`
environment variable.  Do this permanently in your :file:`~/.bashrc`
file::

  export PYTHONPATH=$HOME/ase3k:$PYTHONPATH
  export PATH=$PATH:$HOME/ase3k/tools

or your :file:`~/.cshrc` file::

  setenv PYTHONPATH ${HOME}/ase3k:${PYTHONPATH}
  setenv PATH ${PATH}:${HOME}/ase3k/tools


.. index:: test

Run the tests
=============

Make sure that everything works by running the test suite.  This will
create many files, so run the tests in a new directory::
	
  $ cd
  $ mkdir testase
  $ cd testase
  $ testase.py
  ...
       

If any of the tests fail, then let us know on the `mailing list`_.


.. _mailing list: http://lists.berlios.de/mailman/listinfo/gridpaw-developer


.. index:: License, GPL

License
=======

The CAMPOS Atomic Simulation Environment is released under the GNU
Public License version 2.  See the file LICENSE which accompanies the
downloaded files, or see the license at GNU's web server at
http://www.gnu.org/licenses/gpl.html.
