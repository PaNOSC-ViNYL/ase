import io
import os
import re
import sys
import os.path
import tempfile
import functools

import ase.db
from ase.db.table import Table, all_columns
from ase.visualize import view
from ase.io.png import write_png
from ase.db.summary import Summary

from flask import Flask, render_template, request, send_from_directory


app = Flask(__name__)
connection = None
tables = {}
tmpdir = tempfile.mkdtemp()
next_table_id = 1

# Find numbers in formulas so that we can convert H2O to H<sub>2</sub>O:
SUBSCRIPT = re.compile(r'(\d+)')

                
@app.route('/')
def index():
    global next_table_id
    table_id = int(request.args.get('x', '0'))
    if table_id not in tables:
        table_id = next_table_id
        next_table_id += 1
        query = ''
        columns = list(all_columns)
        sort = 'id'
        limit = 100
        opened = set()
    else:
        query, columns, sort, limit, opened = tables[table_id]

    if 'toggle' in request.args:
        column = request.args['toggle']
        if column in columns:
            columns.remove(column)
            if column == sort.lstrip('-'):
                sort = 'id'
        else:
            columns.append(column)
    elif 'sort' in request.args:
        column = request.args['sort']
        if column == sort:
            sort = '-' + column
        elif '-' + column == sort:
            sort = 'id'
        else:
            sort = column
    elif 'query' in request.args:
        query = request.args['query'].encode()
        limit = int(request.args.get('limit', '0'))
        columns = all_columns
        sort = 'id'
        opened = set()
        
    table = Table(connection)
    table.select(query, columns, sort, limit)
    tables[table_id] = query, table.columns, sort, limit, opened
    table.format(SUBSCRIPT)
    return render_template('table.html', t=table, query=query, sort=sort,
                           limit=limit, tid=table_id, opened=opened)

    
@app.route('/open_row/<int:id>')
def open_row(id):
    table_id = int(request.args['x'])
    opened = tables[table_id][-1]
    if id in opened:
        opened.remove(id)
        return ''
    opened.add(id)
    return render_template('more.html',
                           dct=connection.get(id), id=id, tid=table_id)
    
    
@app.route('/image/<name>')
def image(name):
    path = os.path.join(tmpdir, name).encode()
    if not os.path.isfile(path):
        id = int(name[:-4])
        atoms = connection.get_atoms(id)
        if atoms:
            size = atoms.positions.ptp(0)
            i = size.argmin()
            rotation = ['-90y', '90x', ''][i]
            size[i] = 0.0
            scale = min(20, 20 / size.max() * 10.0)
        else:
            scale = 20
            rotation = ''
        write_png(path, atoms, show_unit_cell=1,
                  rotation=rotation, scale=scale)
    return send_from_directory(tmpdir, name)
    
    
@app.route('/gui/<int:id>')
def gui(id):
    atoms = connection.get_atoms(id)
    view(atoms)
    return '', 204, []
        
        
@app.route('/id/<int:id>')
def summary(id):
    s = Summary(connection.get(id), 'html')
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
                   ]  # ('Content-type', 'application/sqlite3')]
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
    globals()['connection'] = ase.db.connect(sys.argv[1])
    app.run(debug=True)
