import os
import sys
import gettext


def enable_localization():
    """Enables localization using gettext.
    
    Translations will be loaded from mo-files when possible.
    """

    domain = 'ag'
    localedir = '%s/po/' % os.path.dirname(__file__)
    
    gettext.bindtextdomain(domain, localedir)
    gettext.textdomain(domain)
    translation = gettext.translation(domain, localedir, fallback=True)
    if sys.version_info[0] == 2:
        translation.install(unicode=True)
    else:
        translation.install()
