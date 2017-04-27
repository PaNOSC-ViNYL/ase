"""WSGI Flask-app for browsing a database.

You can launch Flask's local webserver like this::

    $ ase db abc.db -w

For a real webserver, you need to set the $ASE_DB_APP_CONFIG environment
variable to point to a configuration file like this::

    ASE_DB_NAMES = ['/path/to/db-file/project1.db',
                    'postgresql://user:pw@localhost:5432/project2']
    ASE_DB_HOMEPAGE = '<a href="https://home.page.dk">HOME</a> ::'

Start with something like::

    twistd web --wsgi=ase.db.app.app --port=8000

"""

from __future__ import print_function
import collections
import functools
import io
import os
import os.path as op
import re
import sys
import tempfile

from flask import Flask, render_template, request, send_from_directory, flash

try:
    import matplotlib
    matplotlib.use('Agg', warn=False)
except ImportError:
    pass

import ase.db
import ase.db.web
from ase.db.plot import atoms2png
from ase.db.summary import Summary
from ase.db.table import Table, all_columns
from ase.visualize import view


default_key_descriptions = {
    'id': ('ID', 'Uniqe row ID', 'int', ''),
    'age': ('Age', 'Time since creation', 'str', ''),
    'formula': ('Formula', 'Chemical formula', 'str', ''),
    'user': ('Username', 'Username', 'str', ''),
    'calculator': ('Calculator', 'ASE-calculator name', 'str', ''),
    'energy': ('Energy', 'Total energy', 'float', 'eV'),
    'fmax': ('Maximum force', 'Maximum force', 'float', 'eV/Ang'),
    'charge': ('Charge', 'Charge', 'float', '|e|'),
    'mass': ('Mass', 'Mass', 'float', 'au'),
    'magmom': ('Magnetic moment', 'Magnetic moment', 'float', 'au'),
    'unique_id': ('Unique ID', 'Unique ID', 'float', ''),
    'volume': ('Volume', 'Volume of unit-cell', 'float', '`Ang^3`')}

# Every client-connetions gets one of these tuples:
Connection = collections.namedtuple(
    'Connection',
    ['project',  # project name
     'query',  # query string
     'nrows',  # number of rows matched
     'page',  # page number
     'columns',  # what columns to show
     'sort',  # what column to sort after
     'limit'])  # number of rows per page

app = Flask(__name__)

app.secret_key = 'asdf'

databases = {}
home = ''  # link to homepage
open_ase_gui = True  # click image to open ASE's GUI

# List of (project-name, title) tuples (will be filled in at run-time):
projects = []


def connect_databases(uris):
    for uri in uris:
        if uri.startswith('postgresql://'):
            project = uri.rsplit('/', 1)[1]
        else:
            project = uri.rsplit('/', 1)[-1].split('.')[0]
        databases[project] = ase.db.connect(uri)


next_con_id = 1
connections = {}

tmpdir = tempfile.mkdtemp()  # used to cache png-files

if 'ASE_DB_APP_CONFIG' in os.environ:
    app.config.from_envvar('ASE_DB_APP_CONFIG')
    connect_databases(app.config['ASE_DB_NAMES'])
    home = app.config['ASE_DB_HOMEPAGE']
    open_ase_gui = False
    try:
        os.unlink('tmpdir')
    except FileNotFoundError:
        pass
    os.symlink(tmpdir, 'tmpdir')

# Find numbers in formulas so that we can convert H2O to H<sub>2</sub>O:
SUBSCRIPT = re.compile(r'(\d+)')


def database():
    return databases[request.args.get('project', 'default')]


def prefix():
    if 'project' in request.args:
        return request.args['project'] + '-'
    return ''


errors = 0


def error(e):
    """Write traceback and other stuff to 00-99.error files."""
    global errors
    import traceback
    x = request.args.get('x', '0')
    try:
        cid = int(x)
    except ValueError:
        cid = 0
    con = connections.get(cid)
    with open(op.join(tmpdir, '{:02}.error'.format(errors % 100)), 'w') as fd:
        print(repr((errors, con, e, request)), file=fd)
        if hasattr(e, '__traceback__'):
            traceback.print_tb(e.__traceback__, file=fd)
    errors += 1
    raise e


app.register_error_handler(Exception, error)


@app.route('/')
def index():
    global next_con_id

    if not projects:
        # First time: initialize list of projects
        projects[:] = [(proj, d.metadata.get('title', proj))
                       for proj, d in sorted(databases.items())]

    con_id = int(request.args.get('x', '0'))
    if con_id in connections:
        project, query, nrows, page, columns, sort, limit = connections[con_id]
        newproject = request.args.get('project')
        if newproject is not None and newproject != project:
            con_id = 0

    if con_id not in connections:
        # Give this connetion a new id:
        con_id = next_con_id
        next_con_id += 1
        project = request.args.get('project', projects[0][0])
        query = ''
        nrows = None
        page = 0
        columns = None
        sort = 'id'
        limit = 25

    db = databases[project]

    if not hasattr(db, 'meta'):
        meta = build_metadata(db)
        db.meta = meta
    else:
        meta = db.meta

    if columns is None:
        columns = meta.get('default_columns') or list(all_columns)

    origquery = query

    if 'sort' in request.args:
        column = request.args['sort']
        if column == sort:
            sort = '-' + column
        elif '-' + column == sort:
            sort = 'id'
        else:
            sort = column
        page = 0
    elif 'query' in request.args:
        query = request.args['query']
        origquery = query
        for special in meta['special_keys']:
            kind, key = special[:2]
            if kind == 'SELECT':
                value = request.args['select_' + key]
                special[-1] = value
                if value:
                    query += ',{}={}'.format(key, value)
            elif kind == 'BOOL':
                value = request.args['bool_' + key]
                special[-1] = value
                if value:
                    query += ',{}={}'.format(key, value)
            else:
                v1 = request.args['from_' + key]
                v2 = request.args['to_' + key]
                var = request.args['range_' + key]
                special[-3:] = [v1, v2, var]
                if v1 or v2:
                    var = request.args['range_' + key]
                    if v1:
                        query += ',{}<={}'.format(v1, var)
                    else:
                        query += ',{}'.format(var)
                    if v2:
                        query += '<={}'.format(v2)
        query = query.lstrip(',')
        sort = 'id'
        page = 0
        nrows = None
    elif 'limit' in request.args:
        limit = int(request.args['limit'])
        page = 0
    elif 'page' in request.args:
        page = int(request.args['page'])

    if 'toggle' in request.args:
        column = request.args['toggle']
        if column == 'reset':
            columns = meta.get('default_columns') or list(all_columns)
        else:
            if column in columns:
                columns.remove(column)
                if column == sort.lstrip('-'):
                    sort = 'id'
                    page = 0
            else:
                columns.append(column)

    okquery = query

    if nrows is None:
        try:
            nrows = db.count(query)
        except (ValueError, KeyError) as e:
            flash(', '.join(['Bad query'] + list(e.args)))
            okquery = 'id=0'  # this will return no rows
            nrows = db.count(okquery)

    table = Table(db)
    table.select(okquery, columns, sort, limit, offset=page * limit)

    con = Connection(project, origquery, nrows, page, columns, sort, limit)
    connections[con_id] = con

    if len(connections) > 1000:
        # Forget old connections:
        for cid in sorted(connections)[:200]:
            del connections[cid]

    table.format(SUBSCRIPT)
    addcolumns = [column for column in all_columns + table.keys
                  if column not in table.columns]

    return render_template('table.html',
                           project=project,
                           projects=projects,
                           t=table,
                           md=meta,
                           con=con,
                           x=con_id,
                           home=home,
                           pages=pages(page, nrows, limit),
                           nrows=nrows,
                           addcolumns=addcolumns,
                           row1=page * limit + 1,
                           row2=min((page + 1) * limit, nrows))


@app.route('/image/<name>')
def image(name):
    id = int(name[:-4])
    name = prefix() + name
    path = op.join(tmpdir, name)
    if not op.isfile(path):
        db = database()
        atoms = db.get_atoms(id)
        atoms2png(atoms, path)

    return send_from_directory(tmpdir, name)


@app.route('/cif/<name>')
def cif(name):
    id = int(name[:-4])
    name = prefix() + name
    path = op.join(tmpdir, name)
    if not op.isfile(path):
        db = database()
        atoms = db.get_atoms(id)
        atoms.write(path)
    return send_from_directory(tmpdir, name)


@app.route('/plot/<png>')
def plot(png):
    png = prefix() + png
    return send_from_directory(tmpdir, png)


@app.route('/gui/<int:id>')
def gui(id):
    if open_ase_gui:
        db = database()
        atoms = db.get_atoms(id)
        view(atoms)
    return '', 204, []


@app.route('/id/<int:id>')
def summary(id):
    db = database()
    prfx = prefix() + str(id)
    s = Summary(db.get(id), db.meta, SUBSCRIPT, prfx, tmpdir)
    return render_template('summary.html',
                           project=request.args.get('project', 'default'),
                           projects=projects,
                           s=s,
                           home=home,
                           md=db.meta,
                           open_ase_gui=open_ase_gui)


def tofile(project, query, type, limit=0):
    fd, name = tempfile.mkstemp(suffix='.' + type)
    con = ase.db.connect(name, use_lock_file=False)
    db = databases[project]
    for row in db.select(query, limit=limit):
        con.write(row,
                  data=row.get('data', {}),
                  **row.get('key_value_pairs', {}))
    os.close(fd)
    data = open(name, 'rb').read()
    os.unlink(name)
    return data


def download(f):
    @functools.wraps(f)
    def ff(*args, **kwargs):
        text, name = f(*args, **kwargs)
        headers = [('Content-Disposition',
                    'attachment; filename="{0}"'.format(name)),
                   ]  # ('Content-type', 'application/sqlite3')]
        return text, 200, headers
    return ff


@app.route('/xyz/<int:id>')
@download
def xyz(id):
    fd = io.StringIO()
    from ase.io.xyz import write_xyz
    db = database()
    write_xyz(fd, db.get_atoms(id))
    data = fd.getvalue()
    return data, '{0}.xyz'.format(id)


@app.route('/json')
@download
def jsonall():
    con_id = int(request.args['x'])
    con = connections[con_id]
    data = tofile(con.project, con.query, 'json', con.limit)
    return data, 'selection.json'


@app.route('/json/<int:id>')
@download
def json1(id):
    project = request.args.get('project', 'default')
    data = tofile(project, id, 'json')
    return data, '{0}.json'.format(id)


@app.route('/sqlite')
@download
def sqliteall():
    con_id = int(request.args['x'])
    con = connections[con_id]
    data = tofile(con.project, con.query, 'db', con.limit)
    return data, 'selection.db'


@app.route('/sqlite/<int:id>')
@download
def sqlite1(id):
    project = request.args.get('project', 'default')
    data = tofile(project, id, 'db')
    return data, '{0}.db'.format(id)


@app.route('/robots.txt')
def robots():
    return 'User-agent: *\nDisallow: /\n', 200


def pages(page, nrows, limit):
    """Helper function for pagination stuff."""
    npages = (nrows + limit - 1) // limit
    p1 = min(5, npages)
    p2 = max(page - 4, p1)
    p3 = min(page + 5, npages)
    p4 = max(npages - 4, p3)
    pgs = list(range(p1))
    if p1 < p2:
        pgs.append(-1)
    pgs += list(range(p2, p3))
    if p3 < p4:
        pgs.append(-1)
    pgs += list(range(p4, npages))
    pages = [(page - 1, 'previous')]
    for p in pgs:
        if p == -1:
            pages.append((-1, '...'))
        elif p == page:
            pages.append((-1, str(p + 1)))
        else:
            pages.append((p, str(p + 1)))
    nxt = min(page + 1, npages - 1)
    if nxt == page:
        nxt = -1
    pages.append((nxt, 'next'))
    return pages


def build_metadata(db):
    meta = db.metadata

    mod = {}
    if db.python:
        with open(db.python) as fd:
            exec(compile(fd.read(), db.python, 'exec'), mod)

    for key, default in [('title', 'ASE database'),
                         ('default_columns', []),
                         ('special_keys', []),
                         ('key_descriptions', {}),
                         ('layout', [])]:
        meta[key] = mod.get(key, meta.get(key, default))

    if not meta['default_columns']:
        meta['default_columns'] = ['id', 'formula']

    # Also fill in default key-descriptions:
    kd = default_key_descriptions.copy()
    kd.update(meta['key_descriptions'])
    meta['key_descriptions'] = kd

    sk = []
    for special in meta['special_keys']:
        kind = special[0]
        if kind == 'SELECT':
            key = special[1]
            choises = sorted({row.get(key) for row in db.select(key)})
            if key in kd:
                longkey = kd[key][0]
            else:
                longkey = key
            special = ['SELECT', key, longkey, choises, choises[0]]
        elif kind == 'BOOL':
            key = special[1]
            if key in kd:
                longkey = kd[key][0]
            else:
                longkey = key
            special = ['BOOL', key, longkey, '']
        else:
            # RANGE
            special = list(special) + ['', '', special[3][0][1]]
        sk.append(special)
    meta['special_keys'] = sk

    if not meta['layout']:
        keys = ['id', 'formula', 'age']
        meta['layout'] = [
            ('Basic properties',
             ['ATOMS', 'CELL'
              ('Key Value Pairs', keys), 'FORCES'])]

    if mod:
        meta['functions'] = ase.db.web.functions[:]
        ase.db.web.functions[:] = []

    sub = re.compile(r'`(.)_(.)`')
    sup = re.compile(r'`(.*)\^(.)`')
    # Convert LaTeX to HTML:
    for key, value in meta['key_descriptions'].items():
        short, long, type, unit = value
        unit = sub.sub(r'\1<sub>\2</sub>', unit)
        unit = sup.sub(r'\1<sup>\2</sup>', unit)
        meta['key_descriptions'][key] = (short, long, type, unit)

    print(meta)
    return meta


if __name__ == '__main__':
    if len(sys.argv) > 1:
        connect_databases(sys.argv[1:])
    app.run(host='0.0.0.0', debug=True)
