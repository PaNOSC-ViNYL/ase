.. _newrelease:

===========
New release
===========

* Make sure all tests pass.

* Go through the git-logs and make sure all important changes since last
  stable release are mentioned in the :ref:`releasenotes`.

* Build the web-page::

      $ cd doc
      $ make clean
      $ make

  and check the generated images with ``make inspect``.

* Update ``__version__`` to ``"x.y.z"`` in :git:`ase/__init__.py`.

* Upload to PyPI::

      $ python3 setup.py sdist
      $ python3 setup.py bdist_wheel
      $ twine upload dist/*

* Push and make a tag "x.y.z".

* Create pull-request for Easy-Build.  The EasyBuild (``.eb``) files
  are in docs/development/easybuild.

  * Rename the old ``.eb`` files, updating the ASE version number.
    There are two ``.eb`` files, one for Python 2.7.12 and one for
    Python 3.5.2.  Use ``git mv`` to rename the files.

  * Edit file file.

    * Update the **version number**.

    * Remove the checksum (delete the entire line in the file).

  * Check the syntax and style, and insert new checksums by running
    these commands::

      eb --check-style ASE-X.Y.X-Python*.eb
      eb --inject-checksums sha256 ASE-X.Y.X-Python*.eb

  * Submit the new files::

      eb --new-pr --pr-commit-message "ASE updated to version X.Y.Z" ASE-X.Y.X-Python*.eb

  * Commit the updated ``*.eb`` files, so they will be part of the
    *next* release.
    
  If the commands above fails, your need to `integrate EasyBuild with github`_.

* Export issues, MR's, ... from GitLab (https://gitlab.com/ase/ase/export)
  and store the tar-file in a safe place.

* Merge *master* into the *web-page* branch (which is used for creating the
  web-page for the stable version).

* Update version numbers in :ref:`news`, :ref:`releasenotes` and
  :ref:`download_and_install` pages.

* Increase the version number and push ("x.y+1.0b1").

* Send announcement email to the ``ase-users`` mailing list.

  Number of commits since last release::

      $ git shortlog -s -n 3.13.0..

.. _`integrate EasyBuild with github`: https://wiki.fysik.dtu.dk/niflheim/EasyBuild_modules#setting-up-github-integration
