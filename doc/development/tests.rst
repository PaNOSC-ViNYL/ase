.. module:: test

================
Testing the code
================

All additions and modifications to ASE should be tested.

.. index:: testase.py

Test scripts should be put in the :svn:`ase/test` directory.  The
scripts in this directory may be run by using the script
:svn:`tools/testase.py` or by using the function:

.. function:: test.test(verbosity=1, dir=None)
    
    Runs the test scripts in :svn:`ase/test`.


.. important::

  When you fix a bug, add a test to the test suite checking that it is
  truly fixed.  Bugs sometimes come back, do not give it a second
  chance!
