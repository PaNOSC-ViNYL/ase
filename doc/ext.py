from ase.utils.sphinx import mol_role, git_role_tmpl, create_png_files


def git_role(role, rawtext, text, lineno, inliner, options={}, content=[]):
    return git_role_tmpl('https://gitlab.com/ase/ase/blob/master/',
                         role,
                         rawtext, text, lineno, inliner, options, content)


def setup(app):
    app.add_role('mol', mol_role)
    app.add_role('git', git_role)
    create_png_files()
