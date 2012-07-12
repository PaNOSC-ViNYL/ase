.. _translate:

Translate ASE
=============

You can contribute by translating the ASE GUI, :file:`ag`, into your language.


How to translate
----------------

If any of the below steps prove difficult, be sure to ask on the
:ref:`developer mailing list <mailing_lists>`.  These steps should work on GNU/Linux.

* Download ASE.
* Go to :file:`ase/gui/po`.  There is a directory of the form :file:`{ll}_{LL}` for each language, where :file:`{ll}` is the `language code`_ and :file:`{LL}` the `country code`_.
* If your language is not there already, run :file:`LANG={ll}_{LL} make init`, substituting the desired language/country codes.
* There should now be a template file for your language :file:`{ll}_{LL}/LC_MESSAGES/ag.po` which can be filled out.

You can edit the po-file with any text editor.  It is easiest with a dedicated po-editor such as gtranslator, poedit or the gettext mode in EMACS (the package ":file:`gettext-el`").  Fill out the missing :file:`msgstr` entries like this::

  #: ../energyforces.py:61
  msgid "Calculate potential energy and the force on all atoms"
  msgstr "Beregn potentiel energi og kræfter på alle atomer"

If a translation in your language already exists, there may be some
"fuzzy" messages.  These are translations that were written
previously, but now need to be reviewed.  Maybe because the original
string has been slightly modified.  Edit as appropriate and remove the
"fuzzy" flag.

There will be a few special constructs such as string substitution
codes :file:`%(number)d` or :file:`%s`.  These should remain unchanged
in the translation as they are replaced by numbers or text at runtime.
An underscore like in :file:`msgid "_File"` indicates that `F` is a
shortcut key.  Conflicting shortcut keys are not a big problem, but
avoid them if you see them.  Finally, some messages may have a lot of
whitespace in them.  This is due to bad programming style; just try to
get approximately the same spacing in your translation.

Check and commit your translation
---------------------------------

* You can check the syntax by running :file:`msgfmt -cv ag.po`.  This will report syntax errors.

* You can test your translation in :file:`ag` by first issuing :file:`make` in :file:`ase/gui/po`, then running :file:`LANG={ll}_{LL}.UTF-8 ag`.

Depending on your operating system, you may need to install
:file:`gettext` or :file:`locales`.

Send the partially or completely translated po-file to the developers
mailing list and ask to have it committed.  In fact, we will be quite thrilled
if you send an e-mail even before you start, and be sure to send one
whenever you have questions.

.. note::

  Certain uncommon languages such as Lojban, Anglo-Saxon or Klingon
  may not be compatible with our current build system.  Please let us
  know if you want to translate ASE into such languages.

Maintaining translations
------------------------

Messages will once in a while be added or changed in the ASE.  Running
:file:`make` in :file:`ase/gui/po` automatically synchronizes all templates with
the messages in the current source tree while maximally reusing the
existing translations.  Some strings may be marked "fuzzy", indicating
that they need review by translators (this happens e.g. if an English
message is changed only slightly).  One can then update the few fuzzy
or untranslated messages.  The obvious time to do this is shortly
before a new stable release.

If you are a committer, please run :file:`make` before committing and
briefly check by running the translated ag that nothing is obviously horrible.

.. _language code: http://www.gnu.org/software/gettext/manual/gettext.html#Language-Codes
.. _country code: http://www.gnu.org/software/gettext/manual/gettext.html#Country-Codes
