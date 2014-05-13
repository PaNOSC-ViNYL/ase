from __future__ import print_function
import io
import os
import sys
import tempfile
import functools
import webbrowser

import ase.db
from ase.io.png import write_png
from ase.db.summary import Summary
from ase.db.table import Table

from flask import Flask, render_template, request, send_from_directory


app = Flask(__name__)
table = None


@app.route('/image/<name>')
def image(name):
    if not os.path.isfile('/tmp/' + name):
        id = int(name[:-4])
        a = table.connection.get_atoms(id)
        write_png('/tmp/' + name, a, show_unit_cell=1)
    return send_from_directory('/tmp', name)
    
    
@app.route('/gui/<int:id>')
def gui(id):
    table.gui(id)
    return '', 204, []
        
        
@app.route('/')
def index():
    if 'moreless' in request.args:
        n = int(request.args['moreless'])
        table.moreless(n)
    elif 'toggle' in request.args:
        key = request.args['toggle']
        table.toggle(key)
    elif 'sort' in request.args:
        column = request.args['sort']
        table.toggle_sort(column)
    elif 'query' in request.args:
        try:
            limit = int(request.args.get('limit'))
        except ValueError:
            limit = None
        table.search(request.args['query'], limit)
        
    table.format('html')

    return render_template('table.html', t=table)

    
@app.route('/summary/<int:id>')
def summary(id):
    s = Summary(table.connection.get(id), 'html')
    return render_template('summary.html', s=s)

    
def tojson(dicts):
    fd = io.BytesIO()
    con = ase.db.connect(fd, 'json', use_lock_file=False)
    writedb(con, dicts)
    return fd.getvalue()
    

def tosqlite(dicts):
    fd, name = tempfile.mkstemp(suffix='.db')
    con = ase.db.connect(name, use_lock_file=False)
    writedb(con, dicts)
    os.close(fd)
    data = open(name).read()
    os.unlink(name)
    return data
    

def writedb(con, dicts):
    for dct in dicts:
        con.write(dct,
                  keywords=dct.get('keywords', []),
                  data=dct.get('data', {}),
                  **dct.get('key_value_pairs', {}))


def download(f):
    @functools.wraps(f)
    def ff(*args, **kwargs):
        text, name = f(*args, **kwargs)
        headers = [('Content-Disposition',
                    'attachment; filename="{0}"'.format(name)),
                   ]  #('Content-type', 'application/sqlite3')]
        return text, 200, headers
    return ff
    
    
@app.route('/json')
@download
def jsonall():
    data = tojson(row.dct for row in table.rows)
    return data, 'selection.json'


@app.route('/json/<int:id>')
@download
def json(id):
    dct = table.connection.get(id)
    data = tojson([dct])
    return data, '{0}.json'.format(id)


@app.route('/sqlite')
@download
def sqliteall():
    data = tosqlite(row.dct for row in table.rows)
    return data, 'selection.db'.format(id)

    
@app.route('/sqlite/<int:id>')
@download
def sqlite(id):
    dct = table.connection.get(id)
    data = tosqlite([dct])
    return data, '{0}.db'.format(id)

    
if __name__ == '__main__':
    con = ase.db.connect(sys.argv[1])
    globals()['table'] = Table(con)
    app.run(debug=True)
