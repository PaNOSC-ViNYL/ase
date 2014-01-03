import os
import operator
from time import time
from random import randint

try:
    from functools import wraps
except ImportError:
    wraps = lambda f: lambda g: g  # PY24

from ase.utils import Lock
from ase.atoms import Atoms
from ase.data import atomic_numbers
from ase.constraints import FixAtoms
from ase.parallel import world, broadcast, DummyMPI
from ase.calculators.calculator import get_calculator, all_properties
from ase.calculators.singlepoint import SinglePointCalculator


T0 = 946681200.0  # January 1. 2000

YEAR = 31557600.0  # 365.25 days

seconds = {'s': 1,
           'm': 60,
           'h': 3600,
           'd': 86400,
           'w': 604800,
           'M': 2629800,
           'y': YEAR}

ops = {'<': operator.lt,
       '<=': operator.le,
       '=': operator.eq,
       '>=': operator.ge,
       '>': operator.gt,
       '!=': operator.ne}


def connect(name, type='extract_from_name', create_indices=True,
            use_lock_file=False):
    if type == 'extract_from_name':
        if name is None:
            type = None
        elif name.startswith('postgres://'):
            type = 'postgres'
        else:
            type = os.path.splitext(name)[1][1:]

    if type is None:
        return NoDatabase()

    if type == 'json':
        from ase.db.json import JSONDatabase
        return JSONDatabase(name, use_lock_file=use_lock_file)
    if type == 'db':
        from ase.db.sqlite import SQLite3Database
        return SQLite3Database(name, create_indices, use_lock_file)
    if type == 'postgres':
        from ase.db.postgresql import PostgreSQLDatabase
        return PostgreSQLDatabase(name, create_indices)
    raise ValueError('Unknown database type: ' + type)


class FancyDict(dict):
    def __getattr__(self, key):
        if key not in self:
            return dict.__getattribute__(self, key)
        value = self[key]
        if isinstance(value, dict):
            return FancyDict(value)
        return value

    def __dir__(self):
        return self.keys()
        

def lock(method):
    @wraps(method)
    def new_method(self, *args, **kwargs):
        if self.lock is None:
            return method(self, *args, **kwargs)
        else:
            #with self.lock: PY24
            #    return method(self, *args, **kwargs)
            self.lock.acquire()
            try:
                return method(self, *args, **kwargs)
            finally:
                self.lock.release()
    return new_method


def parallel(method):
    if world.size == 1:
        return method
        
    @wraps(method)
    def new_method(*args, **kwargs):
        ex = None
        result = None
        if world.rank == 0:
            try:
                result = method(*args, **kwargs)
            except Exception as ex:
                pass
        ex, result = broadcast((ex, result))
        if ex is not None:
            raise ex
        return result
    return new_method


def parallel_generator(generator):
    if world.size == 1:
        return generator
        
    @wraps(generator)
    def new_generator(*args, **kwargs):
        if world.rank == 0:
            for result in generator(*args, **kwargs):
                result = broadcast(result)
                yield result
            broadcast(None)
        else:
            result = broadcast(None)
            while result is not None:
                yield result
    return new_generator


class NoDatabase:
    def __init__(self, filename=None, create_indices=True,
                 use_lock_file=False):
        self.filename = filename
        self.create_indices = create_indices
        if use_lock_file:
            self.lock = Lock(filename + '.lock', world=DummyMPI())
        else:
            self.lock = None

        self.timestamp = None  # timestamp form last write

    @parallel
    @lock
    def write(self, atoms, keywords=[], key_value_pairs={}, data={},
              timestamp=None, **kwargs):
        """Write atoms to database with keywords and key-value pairs.
        
        atoms: Atoms object
            Write atomic numbers, positions, unit cell and boundary
            conditions.  If a calculator is attached, write also already
            calculated properties such as the energy and forces.
        keywords: list of str
            List of keywords.
        key_value_pairs: dict
            Dictionary of key-value pairs.  Values must be strings or numbers.
        data: dict
            Extra stuff (not for searching).
            
        Key-value pairs can also be set using keyword arguments::
            
            connection.write(atoms, name='ABC', frequency=42.0)
            
        """
        
        if timestamp is None:
            timestamp = (time() - T0) / YEAR
        self.timestamp = timestamp
        if atoms is None:
            atoms = Atoms()
        
        kvp = dict(key_value_pairs)  # modify a copy
        kvp.update(kwargs)
        id = self._write(atoms, keywords, kvp, data)
        return id
        
    def _write(self, atoms, keywords, key_value_pairs, data):
        return 1

    @parallel
    @lock
    def reserve(self, *keywords, **key_value_pairs):
        """Write empty row if not already present.
        
        Usage::
            
            id = connection.reserve('keyword1', 'keyword2', ...,
                                    key1=value1, key2=value2, ...)
        
        Write an empty row with the given keywords and key-value pairs and
        return the integer id.  If such a row already exists, don't write
        anything and return None.
        """
        
        for dct in self._select(keywords,
                                [(key, '=', value)
                                 for key, value in key_value_pairs.items()]):
            return None

        self.timestamp = (time() - T0) / YEAR
        id = self._write(Atoms(), keywords, key_value_pairs, {})
        return id
        
    def __delitem__(self, id):
        self.delete([id])
        
    def collect_data(self, atoms):
        dct = atoms2dict(atoms)
        dct['timestamp'] = self.timestamp
        dct['username'] = os.getenv('USER')
        if atoms.calc is not None:
            dct['calculator_name'] = atoms.calc.name.lower()
            dct['calculator_parameters'] = atoms.calc.todict()
            if len(atoms.calc.check_state(atoms)) == 0:
                dct.update(atoms.calc.results)
        return dct

    def get_atoms(self, selection=None, attach_calculator=False,
                  add_additional_information=False, **kwargs):
        """Get Atoms object.
        
        selection: int, str or list
            See the select() method.
        attach_calculator: bool
            Attach calculator object to Atoms object (defaul value is
            False).
        add_additional_information: bool
            Put keywords, key-value pairs and data into Atoms.info
            dictionary.
        """
            
        dct = self.get(selection, fancy=False, **kwargs)
        atoms = dict2atoms(dct, attach_calculator)
        if add_additional_information:
            atoms.info = {}
            for key in ['unique_id', 'keywords', 'key_value_pairs', 'data']:
                if key in dct:
                    atoms.info[key] = dct[key]
        return atoms

    def __getitem__(self, selection):
        if selection == slice(None, None, None):
            return [self[None]]
        return self.get_atoms(selection)

    def get(self, selection=None, fancy=True, **kwargs):
        """Select a single row and return it as a dictionary.
        
        selection: int, str or list
            See the select() method.
        fancy: bool
            return fancy dictionary with keys as attributes (this is the
            default).
        """
        dcts = list(self.select(selection, fancy, limit=2, **kwargs))
        assert len(dcts) == 1
        dct = dcts[0]
        if fancy:
            dct = FancyDict(dct)
        return dct

    @parallel_generator
    def select(self, selection=None, fancy=True, filter=None, explain=False,
               verbosity=1, limit=None, **kwargs):
        """Select rows.
        
        Return iterator with results as dictionaries.  Selection is done
        using key-value pairs, keywords and the special keys:
            
            age, username, calculator, energy, magmom and/or charge.
        
        selection: int, str or list
            Can be:
            
            * an integer id
            * a string like 'key=value', where '=' can also be one of
              '<=', '<', '>', '>=' or '!='.
            * a string like 'keyword'
            * comma separated strings like 'key1<value1,key2=value2,keyword'
            * list of strings or tuples: [('charge', '=', 1)].
        fancy: bool
            return fancy dictionary with keys as attributes (this is the
            default).
        filter: function
            A function that takes as input a dictionary and returns True
            or False.
        explain: bool
            Explain query plan.
        verbosity: int
            Possible values: 0, 1 or 2.
        limit: int or None
            Limit selection.
        """
        
        if selection is None:
            expressions = []
        elif isinstance(selection, int):
            expressions = [('id', '=', selection)]
        elif isinstance(selection, list):
            expressions = selection
        else:
            expressions = selection.split(',')
        keywords = []
        comparisons = []
        for expression in expressions:
            if not isinstance(expression, str):
                comparisons.append(expression)
                continue
            if expression.count('<') == 2:
                value, expression = expression.split('<', 1)
                if expression[0] == '=':
                    op = '>='
                    expression = expression[1:]
                else:
                    op = '>'
                key = expression.split('<', 1)[0]
                comparisons.append((key, op, value))
            for op in ['!=', '<=', '>=', '<', '>', '=']:
                if op in expression:
                    break
            else:
                keywords.append(expression)
                continue
            key, value = expression.split(op)
            comparisons.append((key, op, value))

        cmps = [(key, '=', value) for key, value in kwargs.items()]
        for key, op, value in comparisons:
            if key == 'age':
                key = 'timestamp'
                op = {'<': '>', '<=': '>=', '>=': '<=', '>': '<'}[op]
                value = (time() - T0) / YEAR - time_string_to_float(value)
            elif key in atomic_numbers:
                key = atomic_numbers[key]
                value = int(value)
            else:
                try:
                    value = float(value)
                except ValueError:
                    assert op == '='
            cmps.append((key, op, value))

        for dct in self._select(keywords, cmps, explain=explain,
                                verbosity=verbosity, limit=limit):
            if filter is None or filter(dct):
                if fancy:
                    dct = FancyDict(dct)
                yield dct
                
    @lock
    @parallel
    def update(self, ids, add_keywords=[], **add_key_value_pairs):
        m = 0
        n = 0
        for id in ids:
            dct = self._get_dict(id)
            keywords = dct.get('keywords', [])
            for keyword in add_keywords:
                if keyword not in keywords:
                    keywords.append(keyword)
                    m += 1
            key_value_pairs = dct.get('key_value_pairs', {})
            n -= len(key_value_pairs)
            key_value_pairs.update(add_key_value_pairs)
            n += len(key_value_pairs)
            self.timestamp = dct['timestamp']
            self._write(dct, keywords, key_value_pairs,
                        data=dct.get('data', {}))
        return m, n


def atoms2dict(atoms):
    data = {
        'numbers': atoms.numbers,
        'pbc': atoms.pbc,
        'cell': atoms.cell,
        'positions': atoms.positions,
        'unique_id': '%x' % randint(16**31, 16**32 - 1)}
    if atoms.has('magmoms'):
        data['magmoms'] = atoms.get_initial_magnetic_moments()
    if atoms.has('charges'):
        data['charges'] = atoms.get_initial_charges()
    if atoms.has('masses'):
        data['masses'] = atoms.get_masses()
    if atoms.has('tags'):
        data['tags'] = atoms.get_tags()
    if atoms.has('momenta'):
        data['momenta'] = atoms.get_momenta()
    if atoms.constraints:
        data['constraints'] = [c.todict() for c in atoms.constraints]
    return data


def dict2atoms(dct, attach_calculator=False):
    constraint_dicts = dct.get('constraints')
    if constraint_dicts:
        constraints = []
        for c in constraint_dicts:
            assert c.pop('__name__') == 'ase.constraints.FixAtoms'
            constraints.append(FixAtoms(**c))
    else:
        constraints = None

    atoms = Atoms(dct['numbers'],
                  dct['positions'],
                  cell=dct['cell'],
                  pbc=dct['pbc'],
                  magmoms=dct.get('magmoms'),
                  charges=dct.get('charges'),
                  tags=dct.get('tags'),
                  masses=dct.get('masses'),
                  momenta=dct.get('momenta'),
                  constraint=constraints)

    if attach_calculator:
        atoms.calc = get_calculator(dct['calculator_name'])(
            **dct['calculator_parameters'])
    else:
        results = {}
        for prop in all_properties:
            if prop in dct:
                results[prop] = dct[prop]
        if results:
            atoms.calc = SinglePointCalculator(atoms, **results)

    return atoms


def time_string_to_float(s):
    if isinstance(s, (float, int)):
        return s
    s = s.replace(' ', '')
    if '+' in s:
        return sum(time_string_to_float(x) for x in s.split('+'))
    if s[-2].isalpha() and s[-1] == 's':
        s = s[:-1]
    i = 1
    while s[i].isdigit():
        i += 1
    return seconds[s[i:]] * int(s[:i]) / YEAR


def float_to_time_string(t):
    t *= YEAR
    for s in 'yMwdhms':
        if t / seconds[s] > 5:
            break
    return '%d%s' % (round(t / seconds[s]), s)
