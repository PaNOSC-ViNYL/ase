import apsw, os, pickle, time
import numpy as np
from Scientific.IO.NetCDF import NetCDFFile as netCDF

from pyspglib import spglib
from symmetry.sgroup import SpaceGroup as SGROUP

from Cheetah.Template import Template
from ase.calculators.jacapo import *
import tempfile

printsql = False

DB = 'db.sqlite'

setup_sql='''
CREATE TABLE dacapo (
	id INTEGER PRIMARY KEY,
        calculation_type TEXT, -- dacapo_run, bfgs, etc... maps onto a script to run the job
	path TEXT, --absolute path to file
	date_added FLOAT, --time since epoch
        composition TEXT,
	xc TEXT,
	pw FLOAT,
	dw FLOAT,
	ft FLOAT,
	nbands INTEGER,
	nkpts INTEGER,  -- number of BZkpoints
	spinpol BOOLEAN,
	dipole BOOLEAN,
	symmetry BOOLEAN,
	has_ados BOOLEAN,
	has_stress BOOLEAN,
        constraints TEXT, --pickle string from constraint object
	a FLOAT,
	b FLOAT,
	c FLOAT,
	alpha,
	beta FLOAT,
	gamma FLOAT,
        description TEXT,
	CHECK (symmetry IN (0, 1)),
	CHECK (spinpol IN (0, 1)),
	CHECK (has_ados IN (0, 1)),
	CHECK (dipole IN (0, 1)),
	CHECK (has_stress IN (0, 1)),
        UNIQUE(path)
	);

-- this table holds the positions of each atom and the forces on each atom
CREATE TABLE positions (
	id INTEGER PRIMARY KEY,
	dacapo_id INTEGER references dacapo(id) on delete cascade deferrable initially deferred,
	atom_index INTEGER, -- index in the atoms object
	chemical_symbol TEXT,
        magnetic_moment FLOAT,
        atom_tag INTEGER, --stored in atoms object
	pseudopotential TEXT,
        pseudopotential_md5 TEXT,
	x FLOAT,
	y FLOAT,
	z FLOAT
	);

CREATE TABLE tagged_calculations (
	id INTEGER PRIMARY KEY,
        tag TEXT,
	dacapo_id INTEGER NOT NULL references dacapo(id) on delete cascade deferrable initially deferred
);


'''

if not os.path.exists(DB):
    # we need to initialize a database
    conn = apsw.Connection(DB)
    with conn:
        c = conn.cursor()
        c.execute(setup_sql)

def get_formula_string(atoms):
    '''create formula string from an atoms which the calculation will be filed
    under.

    symbols are alphabetized.'''
    f = {}
    for atom in atoms:
        sym = atom.symbol
        if sym in f:
            f[sym] += 1
        else:
            f[sym] = 1

    keys = f.keys()
    keys.sort()

    formula = ''
    for k in keys:
        formula += '%s%i' % (k,f[k])

    return formula

def insert_database_entry(calc):

    atoms = calc.get_atoms()
    path = calc.get_absolute_path()

    conn = apsw.Connection(DB)

    data = {}

    # remove entries that already exist
    with conn:
        db = conn.cursor()
        db.execute('pragma foreign_keys=on;')
        db.execute('select id from dacapo where path=?',(path,))
        results = db.fetchall()
        if len(results) == 1:
            #one id was found, we will reuse it.
            data['id'] = results[0][0]
        else:
            data['id'] = None

        #now we delete everything related to the uuid and path
        db.execute('delete from dacapo where path=?',(path,))
        if data['id']:
            db.execute('delete from positions where dacapo_id=?',(data['id'],))
            db.execute('delete from tagged_calculations where dacapo_id=?',(data['id'],))

    data['path'] = calc.get_absolute_path()
    data['date_added'] = time.time()
    data['composition'] = get_formula_string(atoms)
    data['xc'] = calc.get_xc()
    data['pw'] = float(calc.get_pw())
    data['dw'] = float(calc.get_dw())
    data['ft'] = float(calc.get_ft())
    data['nbands'] = int(calc.get_nbands())
    data['nkpts'] = len(calc.get_bz_k_points())
    data['spinpol'] = calc.get_spinpol()
    dipole = calc.get_dipole()
    if type(dipole) == type({}):
        if 'status' in dipole:
            data['dipole'] = True
        else:
            data['dipole'] = None
    data['symmetry'] = calc.get_symmetry()

    nc = netCDF(calc.get_nc())
    if 'TotalStress' in nc.variables:
        data['has_stress'] = True
    else:
        data['has_stress'] = None

    if ('PrintAtomProjectedDOS' in nc.variables
        and 'AtomProjectedDOS_OrdinalMap' in nc.variables
        and 'AtomProjectedDOS_EnergyResolvedDOS' in nc.variables):
        data['has_ados'] = True
    else:
        data['has_ados'] = None

    constraints = atoms._get_constraints()
    if constraints != []:
        data['constraints'] = pickle.dumps(constraints)
    else:
        data['constraints'] = None

    data['description'] = calc.get_description()

    cell = atoms.get_cell()

    # get a,b,c,alpha,beta,gamma
    from Scientific.Geometry import Vector
    A,B,C = [Vector(x) for x in cell]

    radian = 1
    degree = 2*np.pi*radian/360

    data['a'] = A.length()
    data['b'] = B.length()
    data['c'] = C.length()
    data['alpha'] = B.angle(C)/degree
    data['beta'] =  A.angle(C)/degree
    data['gamma'] = A.angle(B)/degree

    with conn:
        db = conn.cursor()
        db.execute('pragma foreign_keys=on;')

        db.execute('''
insert into dacapo (id,
                    path,
                    date_added,
                    composition,
                    xc,pw,dw,ft,nbands,nkpts,
                    spinpol,dipole,symmetry,has_ados,has_stress,
                    constraints,
                    a,b,c,alpha,beta,gamma,
                    description)
                    values (
                    :id,
                    :path,
                    :date_added,
                    :composition,
                    :xc,:pw,:dw,:ft,:nbands,:nkpts,
                    :spinpol,:dipole,:symmetry,:has_ados,:has_stress,
                    :constraints,
                    :a,:b,:c,:alpha,:beta,:gamma,
                    :description);''',data)

        db.execute('select id from dacapo where path=?',(data['path'],))
        dacapo_id = db.fetchall()[0][0]

        forces = nc.variables.get('DynamicAtomForces', None)

        for i,atom in enumerate(atoms):
            atom_data = {}
            atom_data['dacapo_id'] = dacapo_id
            atom_data['atom_index'] = i
            atom_data['chemical_symbol'] = atom.symbol.strip()

            if atom.tag is not None:
                atom_data['atom_tag'] = int(atom.tag)
            else:
                atom_data['atom_tag'] = None

            magmom = atom.magmom
            if magmom is not None:
                atom_data['magnetic_moment'] = float(atom.magmom)
            else:
                atom_data['magnetic_moment'] = None

            atom_data['pseudopotential'] = calc.get_psp(sym=atom.symbol)
            atom_data['pseudopotential_md5'] = None
            x,y,z = float(atom.x),float(atom.y),float(atom.z)
            atom_data['x'] = x
            atom_data['y'] = y
            atom_data['z'] = z

            db.execute('''
insert into positions (dacapo_id,
                       atom_index,
                       chemical_symbol,
                       magnetic_moment,
                       atom_tag,
                       pseudopotential,
                       pseudopotential_md5,
                       x,
                       y,
                       z)
values (:dacapo_id,
        :atom_index,
        :chemical_symbol,
        :magnetic_moment,
        :atom_tag,
        :pseudopotential,
        :pseudopotential_md5,
        :x,
        :y,
        :z)''',atom_data)

    #now add in the tags.
    with conn:
        db = conn.cursor()
        tags = calc.get_tags()
        if tags is not None:
            for tag in tags:
                db.execute('insert into tagged_calculations (tag, dacapo_id) values (?,?)',(tag, dacapo_id))

sql_template = '''\
select path from dacapo as d
inner join positions as p
on p.dacapo_id=d.id
where
d.composition='$formula'
AND d.xc='$calc.get('xc')'
AND d.pw=$calc.get('pw')
AND d.dw=$calc.get('dw')
AND d.ft=$calc.get('ft')
AND d.nkpts=$nkpts
-- unit cell parameters
AND d.a BETWEEN $a-$tolerance AND $a+$tolerance
AND d.b BETWEEN $b-$tolerance AND $b+$tolerance
AND d.c BETWEEN $c-$tolerance AND $c+$tolerance
AND d.alpha BETWEEN $alpha-$tolerance AND $alpha+$tolerance
AND d.beta BETWEEN $beta-$tolerance AND $beta+$tolerance
AND d.gamma BETWEEN $gamma-$tolerance AND $gamma+$tolerance
AND (
-- position and type of each atom
#for $i,$atom in enumerate($atoms)
(p.atom_index=$i
AND p.chemical_symbol='$atom.symbol'
AND p.x BETWEEN $atom.x-$tolerance AND $atom.x+$tolerance
AND p.y BETWEEN $atom.y-$tolerance AND $atom.y+$tolerance
AND p.z BETWEEN $atom.z-$tolerance AND $atom.z+$tolerance)
OR
#end for
1=10) --this just closes the template code and does nothing
#if $constraints is None
AND d.constraints is null
#else
AND d.constraints="$constraints"
#end if
group by path having count(path)==$natoms;
'''


def search_database(calc,tolerance=0.001):
    'tolerance is an absolute tolerance'
    atoms = calc.get_atoms()
    cell = atoms.get_cell()
    from Scientific.Geometry import Vector
    A,B,C = [Vector(x) for x in cell]

    radian = 1
    degree = 2*np.pi*radian/360

    a = A.length()
    b = B.length()
    c = C.length()
    alpha = B.angle(C)/degree
    beta =  A.angle(C)/degree
    gamma = A.angle(B)/degree

    constraints = atoms._get_constraints()
    if constraints != []:
        constraints = pickle.dumps(constraints)
    else:
        constraints = None

    try:
        k1,k2,k3 = calc.get('kpts')
        nkpts = k1*k2*k3
    except ValueError:
        nkpts = len(calc.get('kpts'))

    formula = get_formula_string(atoms)
    natoms = len(atoms)

    sql = Template(sql_template,
                   searchList=[locals()])

    ###########
    if printsql:
        print sql.respond()
    ##########

    conn = apsw.Connection(DB)
    with conn:
        c = conn.cursor()
        c.execute(sql.respond())
        return c.fetchall()

def get_path_from_database(calc,dir='DFT',prefix=None):
    '''returns the path of the ncfile holding the results that would
    be in calculator'''
    results = search_database(calc)
    if len(results) == 1:
        path = results[0][0]
        if os.path.exists(path):
            return path
        else:
            #the file may have been deleted
            calc.set_nc(path)
            calc.write()
            return path
    elif len(results) > 1:
        for r in results:
            print r[0]
        raise Exception, 'More than one result found'
    else:
        # not found, so make a path and insert the calculator
        #then write the calculator out so it exists?
        print 'No result found'

        if not os.path.isdir(dir): os.makedirs(dir)
        if prefix is None:
            prefix = get_formula_string(calc.get_atoms())+'-'
        (fd,fname) = tempfile.mkstemp(dir=dir,prefix=prefix,suffix='.nc')
        os.close(fd)
        os.unlink(fname)
        calc.set_nc(fname)
        calc.write()
        insert_database_entry(calc)
        return fname

def update_from_database(self,dir=None,prefix=None):
    '''monkey patch to update a calculator from the database'''

    nc = get_path_from_database(self,dir,prefix)
    self.set_nc(nc)
    self.update_input_parameters()

if __name__ == '__main__':

    from bulk import *

    a1 = get_cubic_perovskite_atoms('La','Ti','O',4.0)

    calc = Jacapo(pw=350,dw=350,
                  kpts=(4,4,4),
                  xc='PW91',
                  ft=0.1,
                  nbands=320,
                  atoms=a1)

    print get_path_from_database(calc)

    a2 = get_cubic_perovskite_atoms('La','V','O',4.0)
    calc2 = Jacapo(pw=350,dw=350,
                  kpts=(4,4,4),
                  xc='PW91',
                  ft=0.1,
                   nbands=320,
                  atoms=a2)
    print get_path_from_database(calc2)
