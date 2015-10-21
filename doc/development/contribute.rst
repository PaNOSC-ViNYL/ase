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

All work on the source code takes place on https://gitlab.com using Git_.

.. _Git: https://git-scm.com/


Proposed git workflow
---------------------

The workflow described here have two elements to them:

1. The guidelines of project git-branches as a whole.
2. The workflow of the individual developer, i.e. the standard operating
   procedure when cloning, pullin, pushing etc.

Only the latter will covered at the moment. When the developer team agree
on the former it will be added here too. One comment though:

In general:
    
* never work in master branch locally or on GitLab.
* Make a new branch for whatever you are doing, and when you are done, push
  it to your own official repos and from there make a merge request from that
  branch to official master.
* Do *not* ever merge the main master branch into yours.

This way Master (official and local) is kept clean and synchronized.

You can learn the git basics in several places, e.g.
http://git-scm.com/book/en/v2/Git-Basics-Recording-Changes-to-the-Repository


The first steps as a developer
------------------------------

* register as a user on https://gitlab.com
* Find and/or change your user-name in your account setting. You will need it.
* follow directions on https://gitlab.com/help/ssh/README for how to generate
  and add your own public ssh-key
* go to https://gitlab.com/ase/ase and fork the project via the blue 'fork'
  button so that you get your own GitLab repository to play with.

You will now have a fork situated at https://gitlab.com/your-user-name/ase

* From here on:
    
  - ``upstream`` refers to git@gitlab.com:ase/ase and refers to present official site of ase project.
  - ``origin`` refers to git@gitlab.com:your-user-name/ase
  
  So ``upstream/master``, for example, refers to master branch (i.e. trunk or
  developement branch) of official ase project.

* clone the repos from ``origin`` to your local machine::
    
      $ git clone git@gitlab.com:your-user-name/ase.git
      $ cd ase
      
* Add and track ``upstream`` in your local git (done only once for each local
  git repos)::

      $ git remote add upstream git@gitlab.com:ase/ase
      
  You can always check the list of branches of remote repositories that you
  follow in your local git repository::
      
      $ git remote -v


Making changes
--------------

Changes and/or additions can be made both directly in GitLab (for small
changes :ref:`see below section<making-small-changes>`) and on a branch
on your local clone of official repos.
Prefered way using command line on local machine is:

1) First step is to clean up master (since you should not work in master
   it's ok).

  * Fetch information about heads and references from "upsteam" and store it in
    local git used in later merges and checkouts::
        
        $ git fetch upstream
        
  * Jump into local master branch and make sure that it's same as "upstream"
    master branch::
        
        $git checkout master
        $git reset --hard upstream/master
        
    If this is first time there would be no need for hard reset, unless some time
    has passed since the cloning. Still better safe than sorry.

2) Next you can do changes and additions.

  * checkout a (new) local branch with a relevant name. I use the commit to
    enhance the file contribute.rst as an example::
        
        $ git checkout -b add-contribute-rst

  * edit/add the file(s)

  * Stage the files to be committed using ``git add``::
      
        $ git add contribute.rst

  * Check your status::
      
        $ git status

  * Commit the staged changes and add commit message::
      
        $ git commit -m "ENH: Add developer workflow guidlines"
        
    Read the :ref:`commit message
    section<writing-the-commit-message>` guidlines for commit message

  * Push commits to your GitLab repository::
      
        $ git push --set-upstream origin add-contribute-rst

  * Go to gitlab.com/your-user-name/ase <http://gitlab.com/your-user-name/ase>
    and click on '## branches' button (where ## is the number of branches on your
    repo)

  * Find the branch 'add-contributing-rst' and click '+ Merge Request'

  * Provide informative title and more verbose description in the
    body of the Merge Request form

  * Click the green 'Submit new merge request' button

  * Wait for feedback from the developer community and address concerns as
    needed by adding more commits to the 'add-contribute-rst' branch on your
    personal repository and then pushing to your gitlab repository.
    
  * Once the developer community is satisfied with your merge request,
    anyone with push access to gitlab.com/ase/ase <http://gitlab.com/ase/ase>
    can merge your merge request and it will now be part of the master branch

  * If this commit was for a single file, say, then go ahead and remove the
    branch locally and on origin. But wait until the merge-request is approved::
    
      git branch -D add-contribute-rst
      git branch -D origin/add-contribute-rst


.. _making-small-changes:

Making small changes
--------------------

Say you want to fix a typo somewhere. GitLab has an editing feature that
can come in handy. Here are the steps to do that there:
    
* go to https://gitlab.com/your-user-name/ase
* click "Files" and find the file you want to change
* click "Edit" and fix the typo
* click "Merge Requests" and add your change from the master branch
    
At this point someone will take a look at your change and merge it to the
official repository if the change looks good.


.. _writing-the-commit-message:

Writing the commit message
--------------------------

Commit messages should be clear and follow a few basic rules.  Example::

   ENH: add functionality X to ase.<submodule>.

   The first line of the commit message starts with a capitalized acronym
   (options listed below) indicating what type of commit this is.  Then a blank
   line, then more text if needed.  Lines shouldn't be longer than 72
   characters.  If the commit is related to a ticket, indicate that with
   "See #3456", "See ticket 3456", "Closes #3456" or similar.

Describing the motivation for a change, the nature of a bug for bug fixes or
some details on what an enhancement does are also good to include in a commit
message.  Messages should be understandable without looking at the code
changes.  A commit message like ``MAINT: fixed another one`` is an example of
what not to do; the reader has to go look for context elsewhere.

Standard acronyms to start the commit message with are:

:API: an (incompatible) API change
:BLD: change related to building ase
:BUG: bug fix
:DEP: deprecate something, or remove a deprecated object
:DEV: development tool or utility
:DOC: documentation
:ENH: enhancement
:MAINT: maintenance commit (refactoring, typos, etc.)
:REV: revert an earlier commit
:STY: style fix (whitespace, PEP8)
:TST: addition or modification of tests
:REL: related to releasing ase


Code review
===========

Before you start working on a Merge Request, *please* read our
:ref:`python_codingstandard`.

Hopefully someone will look at your changes and give you some
feedback.  Maybe everything is fine and things can be merged to the official
repository right away, but there could also be some more work to do like:

* make it compatible with all supported Pythons (see
  :ref:`download_and_install`).
* write more comments
* fix docstrings
* write a test
* add some documentation

This code review loop is not something we have invented to prevent you from
contributing - it should be viewed as an opportunity for you to learn how to
write code that fits into the ASE codebase.
