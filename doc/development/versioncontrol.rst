.. _versioncontrol:

===========================
Using version control (SVN)
===========================

The version control system used in ASE development is subversion. A thorough
subversion manual can be found at http://svnbook.red-bean.com/, here
is a brief overview of the most basic features needed when developing ASE.

* Checking out the code::

    svn co https://svn.fysik.dtu.dk/projects/ase/trunk ase-svn

  This retrieves the code tree from the subversion repository and places it in
  ``ase-svn`` directory

  To test only the ase version in the ``ase-svn`` directory, you need to 
  set your :envvar:`PYTHONPATH` environment variable ::

    export PYTHONPATH=$HOME/ase-svn (bash)

  or ::

    setenv PYTHONPATH ${HOME}/ase-svn (csh or tcsh)

  assuming that you downloaded ase-svn in your home directory (``$HOME``).
  In order to write documentation you also need 
  (in the directory ``ase-svn/doc``) ::

     cd ase-svn
     python setup.py install
     cd doc
     sphinx-build . _build
  
 
* Updating the working copy of the code (in the directory ``ase-svn``)::

    svn update

* Checking the status of the working copy (in the directory ``ase-svn``)::

    svn stat

  The status about the files which are not in version control can be
  surpassed with the ``-q`` flag, and the status with respect to latest
  additions in server can be checked with the ``-u`` flag.

* Committing the changes to the repository

  Before sending the changes in the working copy to the repository, working
  copy should be updated. After that, the changes can be send with::

    svn commit -m "Message to describe the committed changes"

  If the ``-m`` option is omitted, an editor is opened for writing the
  log message.

* Adding files or directories to version control::

    svn add filename

  If ``filename`` is directory, also all the files within the
  directory are added. Note that ``svn commit`` is required before the
  new files are actually included in version control.

