from __future__ import print_function
import os
import sys
import webbrowser
from urllib import unquote
from wsgiref.simple_server import make_server

from ase.io.png import write_png

from jinja2 import Environment, PackageLoader


def D(*a, **k):
    print(*a, file=sys.stderr, **k)


D('GO')


def run(rows=None, summary=None):
    app = App(rows)

    httpd = make_server('', 8000, app)
    print('Serving HTTP on port 8000...')
    webbrowser.open('http://localhost:8000', new=1)
    httpd.serve_forever()

    
class App:
    def __init__(self, rows):
        self.rows = rows
        env = Environment(loader=PackageLoader('ase.db', '.'))
        self.template = env.get_template('table.html')
        
    def __call__(self, env, s):
        q = env['QUERY_STRING']
        path = env['PATH_INFO']
        D(q,path)
        
        if path == '/style.css':
            s('200 OK', [('Content-type', 'text/css')])
            D(__file__)
            with open(__file__.rsplit('/', 1)[0] + '/style.css') as fd:
                return [fd.read()]
            
        if path.endswith('.png'):
            s('200 OK', [('Content-type', 'image/png')])
            name = path[1:]
            if not os.path.isfile(name):
                id = int(name[:-4])
                a = self.rows.connection.get_atoms(id)
                write_png(name, a, show_unit_cell=1)
            return [open(name, 'rb').read()]
            
        s('200 OK', [('Content-type', 'text/html')])
        args = {}
        if q:
            for k, v in (x.split('=') for x in q.split('&')):
                if v.isdigit():
                    v = int(v)
                args[k] = v
        if 'query' in args:
            query = args['query']
            if isinstance(query, str):
                query = unquote(query)
            self.rows.search(query)
        elif args:
            getattr(self.rows, args.pop('m'))(**args)
        self.rows.format()
        txt = self.template.render(q=self.rows)
        return [txt.encode()]
