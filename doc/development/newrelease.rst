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

* Create pull-request for Easy-Build ...

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
