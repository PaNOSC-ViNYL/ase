.. _writing_documentation_ase:

=====================
Writing documentation
=====================

We use the Sphinx_ tool to generate the documentation (both HTML
and PDF_).  The documentation is stored in SVN as text files in the
:trac:`doc` directory using the reStructuredText_ markup language.

.. _reStructuredText: http://docutils.sourceforge.net/rst.html
.. _Sphinx: http://sphinx.pocoo.org
.. _PDF: ../ase-manual.pdf

Setting up development environment with Vagrant
===============================================

If you contribue documentation together with the
code and tests (and you should) here is how you create a development enviroment
on a virtual CentOS7 guest machine using Vagrant_. It takes about about 10 minutes.

1. install Vagrant_ on a host:

   - on Windows 7:

     - open command prompt (see http://windows.microsoft.com/en-us/windows/command-prompt-faq) and
       create the C:\ase-develop directory::

         cd C:\ase-develop

     - download the following software under **C:\ase-develop**:

       * http://the.earth.li/~sgtatham/putty/latest/x86/putty-0.63-installer.exe
       * https://dl.bintray.com/mitchellh/vagrant/vagrant_1.7.2.msi
       * http://download.virtualbox.org/virtualbox/4.3.26/VirtualBox-4.3.26-98988-Win.exe

     - install the downloaded software (you will need to click in order to agree to install the software)::

         C:\ase-develop>vagrant_1.7.2.msi /passive
         C:\ase-develop>putty-0.63-installer.exe /silent
         C:\ase-develop>VirtualBox-4.3.26-98988-Win.exe --silent --extract --path .
         C:\ase-develop>VirtualBox-4.3.26-98988-MultiArch_amd64.msi /passive

     The machine should reboot.

   - on Ubuntu/Debian:

     - install virtualbox dependencies and vagrant::

         $ sudo apt-get install -y vagrant

     - add the **unprivileged** user to the `vboxusers` group::

         $ whoami=`whoami`&& sudo usermod -a -G vboxusers $whoami

   - on Fedora/RHEL/CentOS:

     - on RHEL6/CentOS6 only: install virtualbox and EPEL repositories::

         $ su -c "yum -y install wget"
         $ su -c "wget http://download.virtualbox.org/virtualbox/rpm/rhel/virtualbox.repo -O /etc/yum.repos.d/virtualbox.repo"
         $ su -c "yum -y install http://download.fedoraproject.org/pub/epel/6/x86_64/epel-release-6-8.noarch.rpm"

     - on Fedora only: install virtualbox repository::

         $ su -c "yum -y install wget"
         $ su -c "wget http://download.virtualbox.org/virtualbox/rpm/fedora/virtualbox.repo -O /etc/yum.repos.d/virtualbox.repo"

     - install virtualbox dependencies and vagrant::

         $ su -c "yum -y install kernel kernel-devel kernel-headers dkms"
         $ su -c "yum -y install https://dl.bintray.com/mitchellh/vagrant/vagrant_1.7.2_x86_64.rpm"

     Reboot in order to get the latest kernel.

     - install virtualbox and add the **unprivileged** user to the `vboxusers` group::

         $ su -c "yum -y install VirtualBox-4.3"
         $ whoami=`whoami`&& su -c "usermod -a -G vboxusers $whoami"

   Logout and login in order for the user to get the `vboxusers` group permissions.

2. as the **unprivileged** user, create the virtual guest::

     $ cd /tmp
     $ mkdir ase-develop
     $ cd ase-develop
     $ wget https://svn.fysik.dtu.dk/projects/ase/trunk/doc/development/Vagrantfile
     $ vagrant up
     $ vagrant ssh -- -X

3. after you log in into the virtual guest in order to build ASE documentation do::

     $ cd /vagrant
     $ svn co https://svn.fysik.dtu.dk/projects/ase/trunk ase
     $ cd ase/doc
     $ PATH=/vagrant/ase/tools:$PATH PYTHONPATH=/vagrant/ase make

.. note::

    The ASE checkout is available on the host under **/tmp/ase-develop/ase**
    and on the guest under **/vagrant/ase**.
    If you need root access on the virtual guest do: **sudo su -**.

4. deploy just built documentation into the virtual guest webserver::

     $ su -c "cp -rpf /vagrant/ase/doc/build/html /var/www"
     $ su -c "chown -R apache.apache /var/www/html"
     $ su -c "systemctl reload httpd.service"

5. you can access the documentation webpage from the host with::

     $ firefox http://localhost:8080


.. _Vagrant: https://www.vagrantup.com/


Installing Docutils and Sphinx
==============================

If you prefer configuring Sphinx manually follow along.

The reStructuredText_ parser that Sphinx needs, is part of the Docutils_
project.  So, we need to install docutils and sphinx (version>= 1.1.3).

.. _Docutils: http://docutils.sourceforge.net/


Other requirements
==================

When building the documentation, a number of png-files are generated_.
For that to work, you need the following installed:

* scipy
* matplotlib
* epydoc
* povray (optional)
* dvipng
* pdflatex
* bibtex
* AUCTex
* fontconfig
* convert (ImageMagick)


.. _using_sphinx:

Using Sphinx
============

.. highlight:: bash

First, you should take a look at the documentation for Sphinx_ and
reStructuredText_.

If you don't already have your own copy of the ASE package, then get
the :ref:`latest_development_release` and install it.

Then :command:`cd` to the :file:`doc` directory and build the html-pages::

  $ cd ~/ase/doc
  $ make

This might take a long time the first time you do it.

.. Note::

   Make sure that you build the Sphinx documentation using the
   corresponding ASE version by setting the environment variables
   :envvar:`$PYTHONPATH` and :envvar:`$PATH`.

Make your changes to the ``.rst`` files, run
:command:`make` again, check the results and if things
look ok, commit::

  $ emacs index.rst
  $ make
  $ firefox build/html/index.html
  $ svn ci -m "..." index.rst

To build a pdf-file, you do this::

  $ make latex


Extensions to Sphinx
====================

.. highlight:: rest

We have a couple of extensions to Sphinx:

**:mol:**

   Use ``:mol:`CH_3OH``` to get :mol:`CH_3OH`.

**:svn:**

   A role for creating a link to a file in SVN.  If you write
   ``:svn:`ase/atoms.py```, you
   will get: :svn:`ase/atoms.py`.

**:trac:**

   A role for creating a link to a file in Trac.  If you write
   ``:trac:`ase/atoms.py```, you
   will get: :trac:`ase/atoms.py`.

**:epydoc:**

   A role for creating a link to the API-documentation generated with
   epydoc_.  If you
   write ``:epydoc:`ase.atoms.Atoms```, you will get:
   :epydoc:`ase.atoms.Atoms`.

**:math:**

   This role is for inline LaTeX-style math.  Example:
   ``:math:`\sin(x_n^2)``` gives you :math:`\sin(x_n^2)`.  This role
   is actually the default for ASE's documentation, so you can leave
   out the ``:math:`` part like here: ```\sin(x_n^2)```.


**.. math::**

   Write displayed LaTeX-style math.  Example::

     .. math:: \frac{1}{1+x^2}

   gives you:

   .. math:: \frac{1}{1+x^2}


.. _epydoc:  http://epydoc.sourceforge.net/


.. _generated:

Running Python code to create figures
=====================================

If you want to include a picture in your page, *you should not* check
in the png-file to our SVN repositoy!  Instead, you should check in
the Python script you used to generate the picture (you can also
generate csv-files or pdf-files like this).  The first line of the
script should look like this::

    # creates: fig1.png, fig2.png, table1.csv

Sphinx will run the script and generate the files that you can
then use in your rst-file.  Examples:

* :ref:`eos`.  Source: :trac:`doc/tutorials/eos/eos.py`,
  :trac:`doc/tutorials/eos/eos.rst`
* :ref:`lattice_constant`.  Source: :trac:`doc/tutorials/lattice_constant.py`,
  :trac:`doc/tutorials/lattice_constant.rst`


reStructedText in emacs
=======================

.. highlight:: common-lisp

For people using emacs, the `reStructuredText extension`_ is highly
recommended. The intallation procedure is described in the top of the
file, but for most people, it is enough to place it in your emacs
load-path (typically ``.emacs.d/``) and add the lines::

  (add-to-list 'load-path "~/.emacs.d")
  (require 'rst)

somewhere in your ``.emacs`` file.

To make the mode auto load for relevant file extension, you can write
something like::

  (setq auto-mode-alist
        (append '(("\\.rst$" . rst-mode)
                  ("\\.rest$" . rst-mode)) auto-mode-alist))

In your ``.emacs`` file.

.. _reStructuredText extension: http://docutils.sourceforge.net/tools/editors/emacs/rst.el

Updating Sphinx
===============

Starting a new project with sphinx requires an initial configuration.
This is achieved by running :command:`sphinx-quickstart`.
When updating from a very old sphinx you may consider
generating new configuration files and updating the old files accordingly.

**Note** that the current project is configured up-to-date,
so if you are "simply" writing the documentation
you **must** skip the :command:`sphinx-quickstart` step
and focus on :ref:`using_sphinx`.

Here is how do you setup the GPAW project with sphinx:

 - :command:`cd` to the :file:`doc` directory,

 - run :command:`sphinx-quickstart`
   and answer the questions (example given for GPAW)::

    > Root path for the documentation [.]:

    > Separate source and build directories (y/N) [n]:

    > Name prefix for templates and static dir [.]: _

    > Project name: GPAW
    > Author name(s): 2008, CAMd et al.
  
    > Project version: 0.5
    > Project release [0.5]:

    > Source file suffix [.rst]:

    > Name of your master document (without suffix) [index]: contents

    > autodoc: automatically insert docstrings from modules (y/N) [n]: y
    > doctest: automatically test code snippets in doctest blocks (y/N) [n]:
    > intersphinx: link between Sphinx documentation of different projects (y/N) [n]: y

   This will create :file:`doc/conf.py` and :file:`doc/contents.rst`.
   Both these files need to be edited further
   (:file:`doc/conf.py` may for example include
   options for ``sphinx.ext.pngmath``)

