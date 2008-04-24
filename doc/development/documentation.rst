=====================
Writing documentation
=====================


We use the Sphinx_ tool to generate the ASE documentation (both HTML
and PDF).  The documetation is stored in SVN as text files in the
:svn:`doc` directory using the reStructuredText_ markup language.



.. _reStructuredText: http://docutils.sf.net/rst.html
.. _Sphinx: http://sphinx.pocoo.org


Installing Docutils and Sphinx
==============================

The reStructuredText_ parser the Sphinx needs is part of the Docutils_ project.  So, we need to install docutils and sphinx::

  XXX
  XXX


.. _Docutils: http://docutils.sf.net


Using Sphinx
============

.. highlight:: bash

If you don't already have your own copy of the ASE package, then get
that first::

  $ svn chockout https://svn.fysik.dtu.dk/projects/ase/trunk ase
  $ cd ase

Then :command:`cd` to the :file:`doc` directory and build the html-pages::

  $ cd doc
  $ mkdir _build
  $ sphinx-build . _build

Make your changes to the ``.rst`` files, run the
:command:`sphinx-build` again, check the results and if things look
ok, commit::

  $ emacs index.rst
  $ sphinx-build . .build
  $ firefox .build/index.html
  $ svn ci -m "..." index.rst



Extensions to Sphinx
====================

We have a couple of extensions to Sphinx:

.. role:: svn

   A role for creating a link to a file in SVN.  If you write
   ``:svn:`ase/atoms.py```, you
   will get: :svn:`ase/atoms.py`.  The implementation of this role is
   here: :svn:`doc/ext.py`.

.. role:: epydoc

   A role for creating a link to the API-documentation generated with
   epydoc_.  If you
   write ``:epydoc:`ase.atoms.Atoms```, you will get:
   :epydoc:`ase.atoms.Atoms`.  The implementation of this role is
   here: :svn:`doc/ext.py`.

.. role:: math

   This role is for inline LaTex-style math.  Example:
   ``:math:`\sin(x_n^2)``` gives you :math:`\sin(x_n^2)`.

.. directive:: math

   Write displayed LaTex-style math.  Example::

     .. math::

        \frac{1}{1+x^2}

   gives you:

   .. math::

      \frac{1}{1+x^2}

The implemantation of the math role and directive is here:
:svn:`doc/mathml.py`.

If you add the line ``.. default-role:: math``, then you can leave out
the ``:math:`` part like here: ```\sin(x_n^2)```.


.. _epydoc:  http://epydoc.sf.net

reStructedText in emacs
=======================

For people using emacs, the `reStructuredText extension`_ is highly
recommended. The intallation procedure is described in the top of the
file, but for most people, it is enough to place it in your emacs load
path (typically ``.emacs.d/``) and add the line::

  (require 'rst)

somewhere in your ``.emacs`` file.

To make the mode auto load for relevant file extension, you can write something like::

  (setq auto-mode-alist
        (append '(("\\.rst$" . rst-mode)
                  ("\\.rest$" . rst-mode)) auto-mode-alist))

.. reStructuredText extension_ http://docutils.sourceforge.net/tools/editors/emacs/rst.el

How does it work?
=================

::
 
  <Directory "/var/www/html/ase">
    AllowOverride All
  </Directory>

  AddType application/xhtml+xml .html
