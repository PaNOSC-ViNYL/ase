=================
How to contribute
=================

Discussion of ASE development takes place on the :ref:`ase-developer
<mailing_lists>` mailing list and on the #ase IRC channel on freenode.

We welcome new developers who would like to help work on improving
ASE.  If you would like to contribute, you should first tell us what
you want to work on.  Use the mailing list for that.


GitLab repository
=================

All work on the source code takes place on https://gitlab.com using git_.

.. _git: https://git-scm.com/


The first steps as a developer
------------------------------

* register as a user on https://gitlab.com
* Find and/or change your user-name in your account setting. You will need it.
* follow directions on https://gitlab.com/help/ssh/README for how to generate
  and add your own public ssh-key
* go to https://gitlab.com/ase/ase and fork the project so that you
  get your own GitLab repository to play with.

You will now have a fork situated at https://gitlab.com/your-user-name/ase


Making small changes
--------------------

Say you want to fix a typo somewhere.  Here are the steps to do that:
    
* go to https://gitlab.com/your-user-name/ase
* click "Files" and find the file you want to change
* click "Edit" and fix the typo
* click "Merge Requests" and add your change from the master branch
    
At this point someone will take a look at your change and merge it to the
official repository if the change looks good.


Larger changes
--------------

.. highlight:: bash

For larger changes you will want to work with the files on your local
computer.  Here are the initial steps:
    
* Upload your public SSH-key to GitLab if you haven't done it
* clone the repository to your computer::
    
      $ git clone git@gitlab.com:your-user-name/ase.git

  where ``your-user-name`` is your GitLab user name

It's a good idea to make a local branch for your work::
    
    $ cd ase
    $ git checkout -b myfix  # checkout new branch

You can learn the git basics in several places, e.g.
http://git-scm.com/book/en/v2/Git-Basics-Recording-Changes-to-the-Repository
    
Edit, test, stage and commit your changes::
    
    $ cd ase
    $ idle atoms.py  # use your favorite text editor here
    $ cd ..
    $ python setup.py test
    $ git diff  # print difference
    $ git status  # print what is modified and/or staged and recommended action
    $ git add ase/atoms.py  # stage atoms.py for commit
    $ git status
    $ git commit -m "Fixed bug"  # commit changes into your local branch copy

Push your changes to your repository::
    
    $ git push --set-upstream origin myfix

and create a Merge Request on GitLab at https://gitlab.com/your-user-name/ase.

If you want to start working on something new, you should switch to your
master branch, pull from the official repository and branch again::
    
    $ git checkout master
    $ git remote add official git@gitlab.com:ase/ase.git  # do this once only
    $ git pull official master
    $ git checkout -b newstuff

You can use ``git remote -v`` to see which repositories you have
registered locally.
    
Code review
===========

Before you start working on a Merge Request, *please* read our
:ref:`python_codingstandard`.

Hopefully someone will look at your changes and give you some
feedback.  Maybe everything is fine and things can be merged to the official repository right away, but there could also be some more work to do like:

* make it compatible with all supported Pythons (see
  :ref:`download_and_install`).
* write more comments
* fix docstrings
* write a test
* add some documentation

This code review loop is not something we have invented to prevent you from
contributing - it should be viewed as an opportunity for you to learn how to
write code that fits into the ASE codebase.
