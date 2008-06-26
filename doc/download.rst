.. _download:

========================
Download and install ASE
========================


Requirements for ASE
====================

The following packages are required for basic ASE functionality:

1) Python_.
2) NumPy_.

.. _Python: http://www.python.org
.. _NumPy: http://www.scipy.org/NumPy


It is highly recommended (but not required) to install also
matplotlib_ and pygtk_.  Most likely, some or all of these are already
be installed on your system.


.. _matplotlib: http://matplotlib.sourceforge.net
.. _pygtk: http://www.pygtk.org




Installation
============

.. index:: Installation

.. highlight:: bash

Get the source code from svn::

  $ cd $HOME
  $ svn checkout https://svn.fysik.dtu.dk/projects/ase/trunk ase
  $ cd ase
	
or from this tarfile (python-ase-3.0.0.tar.gz_)::

  $ cd $HOME
  $ tar xtzf python-ase-3.0.0.tar.gz
  $ mv python-ase-3.0.0 ase
  $ cd ase

.. _python-ase-3.0.0.tar.gz: python-ase-3.0.0.tar.gz

If you have root-permissions, you do this::

  $ sudo python setup.py install


No root-permissions
-------------------
   
Put the directory :file:`$HOME/ase` in your :envvar:`PYTHONPATH`
environment variable, and the directory :file:`$HOME/ase/tools` in
your :envvar:`PATH` environment variable.  Do this permanently in
your :file:`~/.bashrc` file::

  export PYTHONPATH=$HOME/ase:$PYTHONPATH
  export PATH=$PATH:$HOME/ase/tools

or your :file:`~/.cshrc` file::

  setenv PYTHONPATH ${HOME}/ase:${PYTHONPATH}
  setenv PATH ${PATH}:${HOME}/ase/tools

Instead of :envvar:`HOME`, you may use any other directory.

  

.. index:: test

Run the tests
=============

Make sure that everything works by running the :mod:`test
suite <test>`.  This will create many files, so run the tests in a new
directory::
	
  $ cd /tmp
  $ testase.py
  ...
       

If any of the tests fail, then let us know on the :ref:`ml`.




.. index:: License, GPL

License
=======

XXX put this in a file!

The CAMPOS Atomic Simulation Environment is released under the GNU
Public License version 2.  See the file :svn:`LICENSE.txt` which
accompanies the downloaded files, or see the license at GNU's web
server at http://www.gnu.org/licenses/gpl.html.
