.. _newrelease:

===========
New release
===========

* Make sure all tests pass.

* Build the web-page and check the generated images with ``make inspect``.

* Update ``__version__`` in :git:`ase/__init__.py`.

* Upload to PyPI::

      $ python3 setup.py sdist
      $ python3 setup.py bdist_wheel
      $ twine upload dist/*

* Push and make a tag.

* Update :ref:`news`, :ref:`releasenotes` and :ref:`download_and_install` pages.

* Increase the version number and push.

* Send announcement email to the ``ase-users`` mailing list.

  Number of commits since last release::

      $ git shortlog -s -n 3.13.0..
