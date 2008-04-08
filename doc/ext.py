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

    ref = 'https://wiki.fysik.dtu.dk/stuff/ase/' + ref
    set_classes(options)
    node = nodes.reference(rawtext, 'Epydoc:' + text,
                           refuri=ref,
                           **options)
    return [node], []

def fix_sidebar():
    print 'fixing sidebars...',
    for docname, name in [('index', 'News'), ('contents', 'Quick Links')]:
        t = open('.build/%s.html' % docname).read()
        i1 = t.find('<div class="container">')
        if i1 != -1:
            print docname,
            i2 = t.index('</div>', i1)
            i3 = t.index('<div class="sidebarwrapper">', i2)
            i4 = t.index('<h3>This Page</h3>', i3)
            t = (t[:i1] + t[i2 + 6:i3] + '<div class="sidebarwrapper">' +
                 ('<h4>%s</h4>' % name) + 
                 t[i1 + 23:i2] + t[i4:])
            open('.build/%s.html' % docname, 'w').write(t)
    print

def setup(app):
    app.add_role('svn', svn_role)
    app.add_role('epydoc', epydoc_role)
    import atexit
    atexit.register(fix_sidebar)
    
