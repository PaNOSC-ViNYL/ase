from __future__ import print_function
import os
import sys
import webbrowser
from urllib import unquote
from wsgiref.simple_server import make_server

from ase.io.png import write_png
from ase.db.summary import Summary

from jinja2 import Environment, PackageLoader


def D(*a, **k):
    print(*a, file=sys.stderr, **k)

    
def run(table):
    app = App(table)
    httpd = make_server('', 8000, app)
    print('Serving HTTP on port 8000...')
    webbrowser.open('http://localhost:8000', new=1)
    httpd.serve_forever()

    
class App:
    def __init__(self, table):
        self.table = table
        env = Environment(loader=PackageLoader('ase.db', '.'))
        self.template1 = env.get_template('table.html')
        self.template2 = env.get_template('summary.html')
        
    def __call__(self, env, start_response):
        q = env['QUERY_STRING']
        path = env['PATH_INFO']
        
        args = {}
        if q:
            for k, v in (x.split('=') for x in q.split('&')):
                if v.isdigit():
                    v = int(v)
                else:
                    v = unquote(v)
                args[k] = v
                
        if path == '/style.css':
            start_response('200 OK', [('Content-type', 'text/css')])
            return self.css()
    
        if path.endswith('.png'):
            start_response('200 OK', [('Content-type', 'image/png')])
            return self.png(path[1:])
            
        if path.endswith('.ico'):
            return
            
        if path != '/':
            start_response('200 OK', [('Content-type', 'text/html')])
            return self.summary(path[1:])
            
        start_response('200 OK', [('Content-type', 'text/html')])
        
        if 'query' in args:
            self.table.search(args['query'])
        elif 'm' in args:
            method = getattr(self.table, args.pop('m'))
            method(**args)
            
        self.table.format('html')
        txt = self.template1.render(t=self.table)
        return [txt.encode()]
        
    def css(self):
        with open(__file__.rsplit('/', 1)[0] + '/style.css') as fd:
            return [fd.read()]
            
    def png(self, name):
        if not os.path.isfile(name):
            id = int(name[:-4])
            a = self.table.connection.get_atoms(id)
            write_png(name, a, show_unit_cell=1)
        return [open(name, 'rb').read()]
            
    def summary(self, name):
        s = Summary(self.table.connection.get(int(name)), 'html')
        txt = self.template2.render(s=s)
        return [txt.encode()]

        #s('200 OK', [('Content-Disposition',
        #              'attachment; filename="myfile.json"')])
        #return ['{}'.encode()]
