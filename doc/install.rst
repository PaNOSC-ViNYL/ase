.. _download_and_install:

============
Installation
============

Requirements
============

* Python_ 2.7, 3.4-3.6
* NumPy_ (base N-dimensional array package)

Optional:

* For extra functionality: SciPy_ (library for scientific computing)
* For :mod:`ase.gui`: Matplotlib_ (2D Plotting)

.. _Python: http://www.python.org/
.. _NumPy: http://docs.scipy.org/doc/numpy/reference/
.. _SciPy: http://docs.scipy.org/doc/scipy/reference/
.. _Matplotlib: http://matplotlib.org/
.. _PyPI: https://pypi.python.org/pypi/ase
.. _PIP: https://pip.pypa.io/en/stable/


Installation using system package managers
==========================================

Linux
-----

Major GNU/Linux distributions (including Debian and Ubuntu derivatives,
Arch, Fedora, Red Hat and CentOS) have a ``python-ase`` package
available that you can install on your system. This will manage
dependencies and make ASE available for all users.

.. note::
   Depending on the distribution, this may not be the latest
   release of ASE.

Max OSX (Homebrew)
------------------

Mac users may be familiar with Homebrew_.  Before installing ASE with pip_ as
described in the next section. Homebrew's ``python`` package provides an up-
to-date version of Python 2.7.x and sets up ``pip`` for you::

  $ brew install python

.. _Homebrew: http://brew.sh


.. index:: pip
.. _pip installation:


Installation using pip
======================

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

Now you should be ready to use ASE, but before you start, please `run the
tests`_ as described below.


.. note::

    If your OS doesn't have ``numpy``, ``scipy`` and ``matplotlib`` packages
    installed, you can install them with::

        $ pip install --upgrade --user numpy scipy matplotlib


.. _here: https://docs.python.org/3/library/site.html#site.USER_BASE


.. _download:

Installation from source
========================

As an alternative to ``pip``, you can also get the source from a tar-file or
from Git.

:Tar-file:

    You can get the source as a `tar-file <http://xkcd.com/1168/>`__ for the
    latest stable release (ase-3.14.0.tar.gz_) or the latest
    development snapshot (`<snapshot.tar.gz>`_).

    Unpack and make a soft link::

        $ tar -xf ase-3.14.0.tar.gz
        $ ln -s ase-3.14.0 ase

    Here is a `list of tarballs <https://pypi.python.org/simple/ase/>`__.

:Git clone:

    Alternatively, you can get the source for the latest stable release from
    https://gitlab.com/ase/ase like this::

        $ git clone -b 3.14.0 https://gitlab.com/ase/ase.git

    or if you want the development version::

        $ git clone https://gitlab.com/ase/ase.git

:Pip:

    install git master directly with pip::

        $ pip install --upgrade git+https://gitlab.com/ase/ase.git@master

    The ``--upgrade`` ensures that you always reinstall even if the version
    number hasn't changed.


Add ``~/ase`` to your :envvar:`PYTHONPATH` environment variable and add
``~/ase/bin`` to :envvar:`PATH` (assuming ``~/ase`` is where your ASE
folder is).  Alternatively, you can install the code with ``python setup.py
install --user`` and add ``~/.local/bin`` to the front of your :envvar:`PATH`
environment variable (if you don't already have that).

Finally, please `run the tests`_.

.. note::

    We also have Git-tags for older stable versions of ASE.
    See the :ref:`releasenotes` for which tags are available.  Also the
    dates of older releases can be found there.


.. _ase-3.14.0.tar.gz: https://pypi.python.org/packages/d5/2c/
    dd2d0fa59dd45a145151f032cda13c60c49fe1886140d6d36bb444760235/
    ase-3.14.0.tar.gz#md5=ab6e7fd75ca7ad58c9fbc73299eaa7c6


Environment variables
=====================

.. envvar:: PATH

    Colon-separated paths where programs can be found.

.. envvar:: PYTHONPATH

    Colon-separated paths where Python modules can be found.

Set these permanently in your :file:`~/.bashrc` file::

    $ export PYTHONPATH=<path-to-ase-package>:$PYTHONPATH
    $ export PATH=<path-to-ase-command-line-tools>:$PATH

or your :file:`~/.cshrc` file::

    $ setenv PYTHONPATH <path-to-ase-package>:${PYTHONPATH}
    $ setenv PATH <path-to-ase-command-line-tools>:${PATH}

.. note::

   If running on Mac OSX: be aware that terminal sessions will
   source :file:`~/.bash_profile` by default and not
   :file:`~/.bashrc`. Either put any ``export`` commands into
   :file:`~/.bash_profile` or source :file:`~/.bashrc` in all Bash
   sessions by adding

   ::

      if [ -f ${HOME}/.bashrc ]; then
      source ${HOME}/.bashrc
      fi

   to your :file:`~/.bash_profile`.


.. index:: test
.. _running tests:
.. _run the tests:

Test your installation
======================

Before running the tests, make sure you have set your :envvar:`PATH`
environment variable correctly as described in the relevant section above.
Run the tests like this::

    $ ase test  # takes 1 min.

and send us the output if there are failing tests.
