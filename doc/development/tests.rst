================
Testing the code
================

All additions and modifications to ASE should be tested.

Test scripts should be put in the ASE/Tests directory. In this directoory there is a script test.py, that will run all scripts present in both the ASE/Tests and ASE/Examples directories.

The module ASE.Utilities.Tests contains usefull tools for writing tests.

Important: When you fix a bug, add a test to the test suite checking that it is truly fixed. Bugs sometimes come back, do not give it a second chance!

:svn:`ase/test`
