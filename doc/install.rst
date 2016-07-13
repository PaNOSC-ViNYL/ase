.. _download_and_install:

============
Installation
============

Requirements
============

* Python_ 2.6-3.5
* NumPy_ (base N-dimensional array package)

Optional:

* For extra functionality: SciPy_ (library for scientific computing)
* For :mod:`ase.gui`: PyGTK_ (GTK+ for Python) and Matplotlib_ (2D Plotting)

.. _Python: http://www.python.org/
.. _NumPy: http://docs.scipy.org/doc/numpy/reference/
.. _SciPy: http://docs.scipy.org/doc/scipy/reference/
.. _Matplotlib: http://matplotlib.org/
.. _pygtk: http://www.pygtk.org/
.. _PyPI: https://pypi.python.org/pypi/ase
.. _PIP: https://pip.pypa.io/en/stable/


Installation using ``pip``
==========================

.. highlight:: bash

The simplest way to install ASE is to use pip_ which will automatically get
the source code from PyPI_::
    
    $ pip install --upgrade --user ase
    
This will install ASE in a local folder where Python can
automatically find it (``~/.local`` on Unix, see here_ for details).  Some
:ref:`cli` will be installed in the following location:
    
=================  ============================
Unix and Mac OS X  ``~/.local/bin``
Homebrew           ``~/Library/Python/X.Y/bin``
Windows            ``%APPDATA%/Python/Scripts``
=================  ============================

Make sure you have that path in your :envvar:`PATH` environment variable.
    
Now you should be ready to use ASE, but before you start, please run the
tests as described below.

.. note::

    If your OS doesn't have ``numpy``, ``scipy`` and ``matplotlib`` packages
    installed, you can install them with::
        
        $ pip install --upgrade --user numpy scipy matplotlib
    
.. note::

    Some Linux distributions have an ASE package (named ``python-ase``),
    that you can install on your system so that it is avilable for all
    users.

.. _here: https://docs.python.org/3/library/site.html#site.USER_BASE


.. index:: test
.. _running tests:

Test your installation
======================

Run the tests like this::
    
    $ python -m ase.test  # takes 1 min.

and send us the output if there are failing tests.


.. _download:

Installation from source
========================

As an alternative to ``pip``, you can also get the source from a tar-file or
from Git.

:Tar-file:

    You can get the source as a `tar-file <http://xkcd.com/1168/>`__ for the
    latest stable release (ase-3.11.0.tar.gz_) or the latest
    development snapshot (`<snapshot.tar.gz>`_).

    Unpack and make a soft link::
    
        $ tar -xf ase-3.11.0.tar.gz
        $ ln -s ase-3.11.0 ase

:Git clone:

    Alternatively, you can get the source for the latest stable release from
    https://gitlab.com/ase/ase like this::
    
        $ git clone -b 3.11.0 https://gitlab.com/ase/ase.git

    or if you want the development version::

        $ git clone https://gitlab.com/ase/ase.git
    
Add ``~/ase`` to your :envvar:`PYTHONPATH` environment variable and add
``~/ase/tools`` to :envvar:`PATH` (assuming ``~/ase`` is where your ASE
folder is).  Alternatively, you can install the code with ``python setup.py
install --user`` and add ``~/.local/bin`` to the front of your :envvar:`PATH`
environment variable (if you don't already have that).
    
.. note::
    
    We also have Git-tags for older stable versions of ASE.
    See the :ref:`releasenotes` for which tags are available.  Also the
    dates of older releases can be found there.

.. _ase-3.11.0.tar.gz: https://pypi.python.org/packages/fc/7b/
    558e7321f7a879c034ead5d10789b9d6f41beabaee0b156e807c19422ad0/
    ase-3.11.0.tar.gz

    
Environment variables
=====================

.. envvar:: PATH

    Colon-separated paths where programs can be found.
    
.. envvar:: PYTHONPATH

    Colon-separated paths where Python modules can be found.

Set these permanently in your :file:`~/.bashrc` file::

    $ export PYTHONPATH=~/ase:$PYTHONPATH
    $ export PATH=~ase/tools:$PATH

or your :file:`~/.cshrc` file::

    $ setenv PYTHONPATH ${HOME}/ase:${PYTHONPATH}
    $ setenv PATH ${HOME}/ase/tools:${PATH}


Old video tutorial
==================

In the video: Introduction to ASE, followed by a guide to installing ASE on a
Linux system.

.. note::

   Use "Right Click -> Play" to play.

.. raw:: html

        <p></p>
        <object width="800" height="600">
        <embed src="https://wiki.fysik.dtu.dk/ase-files/oi_en_800x600.swf"
        type="application/x-shockwave-flash"
        allowFullScreen="false"
        allowscriptaccess="never"
        loop="false"
        play="false"
        width="800" height="600">
        <p></p>
        Video not playing? Download avi <a href="https://wiki.fysik.dtu.dk/ase-files/oi_en.avi">file</a> instead.
        </embed></object>
        <p></p>
