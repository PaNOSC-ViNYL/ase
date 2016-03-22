from ase.utils.sphinx import mol_role
from ase.utils.sphinx import git_role_tmpl, epydoc_role_tmpl
from ase.utils.sphinx import create_png_files


def git_role(role, rawtext, text, lineno, inliner, options={}, content=[]):
    return git_role_tmpl('https://gitlab.com/ase/ase/blob/master/',
                         role,
                         rawtext, text, lineno, inliner, options, content)


def epydoc_role(role, rawtext, text, lineno, inliner, options={}, content=[]):
    return epydoc_role_tmpl('ase', 'http://wiki.fysik.dtu.dk/ase/epydoc/',
                            role,
                            rawtext, text, lineno, inliner, options, content)


def setup(app):
    app.add_role('mol', mol_role)
    app.add_role('git', git_role)
    app.add_role('epydoc', epydoc_role)
    create_png_files()
