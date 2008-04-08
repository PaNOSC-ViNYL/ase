# -*- coding: utf-8 -*-

from docutils import nodes, utils

from docutils.parsers.rst.roles import set_classes

def svn_role(role, rawtext, text, lineno, inliner, options={}, content=[]):
    ref = 'https://trac.fysik.dtu.dk/projects/ase/browser/trunk/' + text
    set_classes(options)
    node = nodes.reference(rawtext, text, refuri=ref,
                           **options)
    return [node], []

def epydoc_role(role, rawtext, text, lineno, inliner, options={}, content=[]):
    components = text.split('.')
    if components[0] != 'ase':
        components.insert(0, 'ase')
    try:
        for n in range(2, len(components) + 1):
            module = __import__('.'.join(components[:n]))
    except ImportError:
        for component in components[1:n]:
            module = getattr(module, component)
            ref = '.'.join(components[:n])
            if isinstance(module, type):
                ref += '-class.html'
            else:
                ref += '-module.html'
        if n < len(components):
            ref += '#' + components[-1]
    else:
        ref = '.'.join(components) + '-module.html'

    ref = 'https://web2.fysik.dtu.dk/ase/epydoc/' + ref
    set_classes(options)
    node = nodes.reference(rawtext, 'Epydoc:' + text,
                           refuri=ref,
                           **options)
    return [node], []

def setup(app):
    app.add_role('svn', svn_role)
    app.add_role('epydoc', epydoc_role)
    #import atexit
    #atexit.register(fix_sidebar)

"""
def fix_sidebar():
    print 'fixing sidebars...',
    for docname in ['index', 'contents']:
        t = open('.build/%s.html' % docname).read()
        i1 = t.find('<p class="first sidebar-title">')
        if i1 != -1:
            print docname,
            i0 = t.find('<div class="sidebar">')
            i1b = t.find('</p>', i1)
            name = t[i1 + 31:i1b]
            i2 = t.index('</div>', i1)
            i3 = t.index('<div class="sidebarwrapper">', i2)
            i4 = t.index('<h3>This Page</h3>', i3)
            t = (t[:i0] + t[i2 + 6:i3] + '<div class="sidebarwrapper">' +
                 ('<h4>%s</h4>' % name) + 
                 t[i1b + 5:i2] + t[i4:])
            open('.build/%s.html' % docname, 'w').write(t)
    print
"""
    
