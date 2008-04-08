Development
===========


If you would like to contribute to the ASE project, you are most
welcome to do so.  There are many thing on the `To do list`.  Just join the
`developer mailing list` and read below how to code, test and
document your work.


**Coding**
             To maintain a mostly consistent style through all the ASE
             source code, we try to use a common `Code standard`.

**Testing**
             All additions and modifications to ASE should be `tested`.  

**Documenting**

             We suggest that you go easy on the docstring documentation
             in the code and save your strength for the documentation
             to be put in the `ASE Manual`.  Here are some `Documentation guidelines` for
             how to write documentation and where to put it.

**Releasing**          
             There are a number of things to do when `Releasing a new
	     version of ASE`.


Would you like to add an ASE interface to your favorite
calculator?  Here is `how you can do it`.


---------------
Tips and tricks
---------------

readline
--------

Be sure to have these lines in your personal ``.pythonrc`` file::

  import rlcompleter
  import readline
  readline.parse_and_bind("tab: complete")

and point the ``PYTHONSTARTUP`` environment variable at it (see
here_ for details).

.. _here: http://www.python.org/doc/current/lib/module-rlcompleter.html

