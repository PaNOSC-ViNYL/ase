.. _download_and_install:

============
Installation
============

.. contents::
    
Requirements
============

* Python_ 2.6-3.5
* NumPy_ (base N-dimensional array package)

Optional:

* For extra functionality: SciPy_ (library for scientific computing)
* For :mod:`ase.gui`: PyGTK_ (GTK+ for Python) and Matplotlib_ (2D Plotting)

Installation using ``pip``
==========================

.. highlight:: bash

The simplest way to install ASE is::
    
    $ pip install ase
    
This will install ASE in your ``~/.local`` folder where Python can
automatically find it.  Make sure you have ``~/.local/bin`` in your
``$PATH`` environment variable.


.. _running_tests:

Testing
=======

Please run the tests::
    
    $ python -m ase.test  # takes 1 min.

and send us the output if there are failing tests.


Installation
------------

Add ``~/ase`` to your $PYTHONPATH environment variable and add
``~/ase/tools`` to $PATH (assuming ``~/ase`` is where your ASE folder is).
    
    


    
    






For the installation on personal laptops we recommend
the binary packages provided for popular Linux distributions
(:ref:`installationguide_package`)
and MS Windows (:ref:`installationguide_windows`).

Please skip to :ref:`installationguide_manual` if you prefer
to install from sources.

If you are on Mac OSX, please follow :ref:`installationguide_macosx`.


.. _installationguide_package:

Installation with package manager on Linux
==========================================

Install the binaries with the software package manager of your Linux
distribution.

This is **the preferred** way to install on a Linux system.

If you prefer to install from sources follow :ref:`installationguide_manual`.

The currently supported systems include (issue the commands below **as root**):

- Fedora::

    yum install -y python-ase

- RHEL/CentOS - available after enabling https://fedoraproject.org/wiki/EPEL::

    yum install -y python-ase

- openSUSE 13.1::

    zypper ar -f http://download.opensuse.org/repositories/home:/dtufys/openSUSE_13.1/home:dtufys.repo
    yast -i python-ase
    yast -i python-matplotlib # optionally

- Debian 7.0::

    sudo bash -c 'echo "deb http://download.opensuse.org/repositories/home:/dtufys/Debian_7.0 /" > /etc/apt/sources.list.d/home_dtufys.sources.list'
    wget http://download.opensuse.org/repositories/home:/dtufys/Debian_7.0/Release.key && sudo apt-key add Release.key && rm Release.key
    sudo apt-get update
    sudo apt-get -y install python-ase
    sudo apt-get -y install python-matplotlib # optionally

- Ubuntu 14.04::

    sudo bash -c 'echo "deb http://download.opensuse.org/repositories/home:/dtufys/xUbuntu_14.04 /" > /etc/apt/sources.list.d/home_dtufys.sources.list'
    wget http://download.opensuse.org/repositories/home:/dtufys/xUbuntu_14.04/Release.key && sudo apt-key add Release.key && rm Release.key
    sudo apt-get update
    sudo apt-get -y install python-ase
    sudo apt-get -y install python-matplotlib # optionally

For the full list of supported distributions check
https://build.opensuse.org/package/show?package=python-ase&project=home%3Adtufys

After performing the installation do not forget to :ref:`running_tests`!


.. _installationguide_macosx:

Installation on OS X
====================

For installation with http://brew.sh please follow
instructions at the `Homebrew ASE installation page
<https://wiki.fysik.dtu.dk/gpaw/install/MacOSX/homebrew.html>`_.

After performing the installation do not forget to :ref:`running_tests`!




.. _installationguide_manual:

Installation from source
========================

ASE binaries are available only for the :ref:`latest_stable_release`,
and all available ASE releases are listed at the :ref:`download` page.

If you need a development version (or a historic version) of ASE
perform a manual installation according to instructions below.
Follow the same instructions if you are configuring ASE on an HPC cluster.

This is the **preferred** way of manually installing ASE.
It offers the following advantages:

- installation is limited to standard user's account:
  it does not pollute the root filesystem,

- user gains access to version control updates, if necessary.



Installation process
--------------------

After the :ref:`download` of ASE source create the link
to the requested version, e.g.:

- if retrieved from SVN::

   $ cd ~
   $ ln -s ase-3.9.1 ase
    
- if retrieved as tar-file::

   $ cd ~
   $ tar -xf python-ase-3.9.1.4567.tar.gz
   $ ln -s python-ase-3.9.1.4567 ase

It is sufficient to
put the directory :file:`$HOME/ase` in your :envvar:`PYTHONPATH`
environment variable, and the directory :file:`$HOME/ase/tools` in
your :envvar:`PATH` environment variable.  Do this permanently in
your :file:`~/.bashrc` file::

  export PYTHONPATH=$HOME/ase:$PYTHONPATH
  export PATH=$HOME/ase/tools:$PATH

or your :file:`~/.cshrc` file::

  setenv PYTHONPATH ${HOME}/ase:${PYTHONPATH}
  setenv PATH ${HOME}/ase/tools:${PATH}

Instead of :envvar:`HOME`, you may use any other directory.

Alternatively, you can install ASE to the user-specific site-packages
directory with::

  $ cd ase
  $ python setup.py install --user

This way, the ASE modules are found on the python path without any
explicit configuration, though you still need to ensure that
:file:`$HOME/.local/bin` (or on Windows,
:file:`%APPDATA%/Python/Scripts`) is on your :envvar:`PATH`.

.. index:: test


.. _download:

Download
--------

.. highlight:: bash

.. _latest_stable_release:

Latest stable release
+++++++++++++++++++++

The latest stable release can be obtained from SVN or as a
`tar-file <http://xkcd.com/1168/>`__.

.. note::

   The recommended installation path is :envvar:`$HOME`.

When using svn please set the following variable:

- bash::

   export ASE_TAGS=https://svn.fysik.dtu.dk/projects/ase/tags/

- csh/tcsh::

   setenv ASE_TAGS https://svn.fysik.dtu.dk/projects/ase/tags/

=======  ===========
Release  Date       
=======  ===========
  3.9.1  Jul 21 2015
  3.9.0  May 28 2015
  3.8.1  Nov 22 2013
  3.8.0  Oct 22 2013
  3.7.1  May 16 2013
  3.7.0  May 13 2013
  3.6.0  Feb 24 2012
  3.5.1  May 24 2011
  3.4.1  Aug 11 2010
  3.4.0  Apr 23 2010
  3.3.1  Jan 20 2010
  3.2.0  Sep  4 2009
  3.1.0  Mar 27 2009
  3.0.0  Nov 13 2008
=======  ===========



.. _latest_development_release:

Latest development release
++++++++++++++++++++++++++

The latest revision can be obtained like this::

  $ git clone https://gitlab.com/ase/ase.git

or from the daily snapshot: `<snapshot.tar.gz>`_.


Video tutorial
==============

In the video: :ref:`overview` of the features of ASE,
followed by a :ref:`installationguide_manual` of ASE on a Linux system.

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

        
.. _installationguide_windows:

Installation on Windows
=======================

.. note::

   ASE is not yet fully functional on Windows!
   https://trac.fysik.dtu.dk/projects/ase/ticket/62

Python(x,y), on both 32- and 64-bit Windows,
requires Microsoft Visual C++ 2008 Redistributable Package (x86),
download and install it from:
https://www.microsoft.com/en-us/download/details.aspx?id=5582
Use http://www.dependencywalker.com/ to find missing DLLs in case of
"ImportError: DLL load failed: The specified module could not be found".

Continue with:

.. note:: installation assumes the python TARGETDIR C:\\Python27,
          leave also the default C:\\Program Files\\pythonxy.

-  pythonxy_. Download the *2.7.5.2* exe installer (other versions
   may be incompatible)and install with::

     Python(x,y)-2.7.5.2.exe /Log="%TMP%\pythonxy_install.log" /S

.. note::

   Open Task Manager and control when the process in finished.

- pygtk_win32_. Download the msi **pygtk-all-in-one** installer.
  Specify the correct TARGETDIR and install::

     pygtk-all-in-one-2.24.2.win32-py2.7.msi TARGETDIR="%HOMEDRIVE%\Python27" ALLUSERS=1 /l*vx "%TMP%\pygtk_install.log" /passive

.. note::

   If performing clicking-installation make sure that the default
   python Windows TARGETDIR is selected.

- Download the python-ase-win32.msi_ installer and install with::

     python-ase-X.X.X.win32.msi /l*vx "%TMP%\python-ase_install.log" /passive

.. note::

   You can build the msi ASE package on Windows with::

      python setup.py bdist_msi

   The msi package will be created under the *dist* directory.

.. _pythonxy: http://code.google.com/p/pythonxy
.. _pygtk_win32: http://ftp.gnome.org/pub/GNOME/binaries/win32/pygtk/2.24/

.. _python-ase-win32.msi:
    https://wiki.fysik.dtu.dk/ase-files/python-ase.win32.msi

After performing the installation do not forget to :ref:`running_tests`!


.. _Python: http://www.python.org/
.. _NumPy: http://docs.scipy.org/doc/numpy/reference/
.. _SciPy: http://docs.scipy.org/doc/scipy/reference/
.. _Matplotlib: http://matplotlib.org/
.. _pygtk: http://www.pygtk.org/

.. _python-ase-3.9.1.4567.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.9.1.4567.tar.gz
.. _python-ase-3.9.0.4465.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.9.0.4465.tar.gz
.. _python-ase-3.8.1.3440.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.8.1.3440.tar.gz
.. _python-ase-3.8.0.3420.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.8.0.3420.tar.gz
.. _python-ase-3.7.1.3184.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.7.1.3184.tar.gz
.. _python-ase-3.7.0.3168.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.7.0.3168.tar.gz
.. _python-ase-3.6.0.2515.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.6.0.2515.tar.gz
.. _python-ase-3.5.1.2175.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.5.1.2175.tar.gz
.. _python-ase-3.4.1.1765.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.4.1.1765.tar.gz
.. _python-ase-3.4.0.1574.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.4.0.1574.tar.gz
.. _python-ase-3.3.1.1390.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.3.1.1390.tar.gz
.. _python-ase-3.2.0.1121.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.2.0.1121.tar.gz
.. _python-ase-3.1.0.846.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.1.0.846.tar.gz
.. _python-ase-3.0.0.657.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.0.0.657.tar.gz
