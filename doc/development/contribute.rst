=================
How to contribute
=================

Discussion of ASE development takes place on the :ref:`ase-developer
<mailing_lists>` mailing list, and also sometimes on the #gpaw IRC
channel on freenode.

We welcome new developers who would like to help work on improving
ASE.  If you would like to contribute, your should first tell us what
you want to work on.  Use the mailing list for that.


SVN access
==========

We don't give new contributers write access to our SVN repository from
day one.  So, you will have to create a patch and send it to the
mailing list::

  $ svn checkout https://svn.fysik.dtu.dk/projects/ase/trunk myase
  $ cd myase
  $ # do your thing ...
  $ svn diff > patch.txt

One of the current committers will look at the patch and give you
some feedback.  Maybe the patch is fine and the committer will commit
it to trunk.  There could also be some more work to do like:

* write more comments
* follow our :ref:`python_codingstandard`
* fix this docstring
* write a test
* add some documentation

Once everyone is happy, the patch can be appied.  This patch-feedback
loop is not something we have invented to prevent you from
contributing - it should be viewed as an opportunity for you to learn
how to write good code that fits into the ASE codebase.  

After a couple of contributions, we will probably trust you enough to
add you as a committer.


Committers
==========

Here is the list of current committers:

==========  ====================================
user name   real name
==========  ====================================
anro        Anthony Goodrow
askhl       Ask Hjorth Larsen
bjork       Jonas Bjork
carstenr    Carsten Rostgaard
dlandis     David Landis
dulak       Marcin Dulak
getri       George Tritsaris
grabow      Lars Grabow
hahansen    Heine Anton Hansen
hanke       Felix Hanke
jakobb      Jakob Blomquist
jblomqvist  Janne Blomqvist
jensj       Jens Jørgen Mortensen
jingzhe     Jingzhe Chen
jkitchin    John Kitchin
jussie      Jussi Enkovaara
kwj         Karsten Wedel Jacobsen
markus      Markus Kaukonen
miwalter    Michael Walter
moses       Poul Georg Moses
mvanin      Marco Vanin
s032082     Christian Glinsvad
s042606     Janosch Michael Rauba
s052580     Troels Kofoed Jacobsen
schiotz     Jakob Schiotz
strange     Mikkel Strange
tjiang      Tao Jiang
==========  ====================================

.. epost={'askhl': 'askhl fysik,dtu,dk', 'mvanin': 'mvanin fysik,dtu,dk', 'kwj': 'kwj fysik,dtu,dk', 'grabow': 'grabow fysik,dtu,dk', 'jussie': 'jussi,enkovaara csc,fi', 'dulak': 'dulak fysik,dtu,dk', 'anro': 'anro fysik,dtu,dk', 'bjork': 'J,Bjork liverpool,ac,uk', 'jkitchin': 'jkitchin andrew,cmu,edu', 'jensj': 'jensj fysik,dtu,dk', 'hanke': 'F,Hanke liverpool,ac,uk', 'jingzhe': 'jingzhe fysik,dtu,dk', 'strange': 'strange fysik,dtu,dk', 'schiotz': 'schiotz fysik,dtu,dk', 'moses': 'poulgeorgmoses gmail,com', 's032082': 's032082 fysik,dtu,dk', 'jblomqvist': 'Janne,Blomqvist tkk,fi', 'dlandis': 'dlandis fysik,dtu,dk', 'getri': 'getri fysik,dtu,dk', 's052580': 's052580 fysik,dtu,dk', 's042606': 's042606 fysik,dtu,dk', 'hahansen': 'hahansen fysik,dtu,dk', 'miwalter': 'Michael,Walter fmf,uni-freiburg,de', 'markus': 'markus,kaukonen iki,fi', 'jakobb': 'jakobb fysik,dtu,dk', 'carstenr': 'carstenr fysik,dtu,dk', 'tjiang': 'tjiang fysik,dtu,dk'}
