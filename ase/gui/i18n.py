import os
import sys
import locale
import gettext

# This module enables localization using gettext.  If this module has
# been imported, translations will be loaded from mo-files when
# possible.

domain = 'ag'
localedir = '%s/po/' % os.path.dirname(__file__)

gettext.bindtextdomain(domain, localedir)
gettext.textdomain(domain)
translation = gettext.translation(domain, localedir, fallback=True)
if sys.version_info[0] == 2:
    translation.install(unicode=True)
else:
    translation.install()
