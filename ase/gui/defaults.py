""" This is a module to handle generic ASE (gui) defaults from a ~/.ASE.rc configuration file, if it exists.
It is imported when opening ag and can then be modified at runtime, if necessary.
syntax for each entry:

key = options
"""

gui_default_settings = {
    'gui_graphs_string' : 'i, e - E[-1]',   # default for the graph command in the gui
    }

def read_defaults():
    import os.path
    name = os.path.expanduser('~/.ASE.rc')
    if os.path.exists(name):
        fd = open(name)
        lines = fd.readlines()
        fd.close()
        for line in lines:
            keys = line.split()
            if gui_default_settings.has_key(keys[0]):
                s = ''
                for k in keys[2:]:
                    s += k + ' '
                gui_default_settings[keys[0]] = s

