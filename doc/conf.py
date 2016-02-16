import sys
import sphinx_rtd_theme

sys.path.append('.')
assert sys.version >= '2.7'

extensions = ['ext',
              'images',
              'sphinx.ext.autodoc',
              'sphinx.ext.mathjax',
              'sphinx.ext.viewcode',
              'sphinx.ext.intersphinx']
source_suffix = '.rst'
#master_doc = 'contents'
project = 'ASE'
copyright = '2016, ASE-developers'
#templates_path = ['templates']
exclude_patterns = ['build']
default_role = 'math'
pygments_style = 'sphinx'
autoclass_content = 'both'
html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
#html_style = 'ase.css'
html_logo = 'static/ase.ico'
html_favicon = 'static/ase.ico'
html_static_path = ['static']
html_last_updated_fmt = '%a, %d %b %Y %H:%M:%S'
intersphinx_mapping = {'gpaw': ('http://wiki.fysik.dtu.dk/gpaw', None),
                       'python': ('http://docs.python.org/2.7', None)}
