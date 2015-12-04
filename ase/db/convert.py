import optparse
import os

from ase.db import connect
from ase.db.sqlite import index_statements


def convert(name):
    con1 = connect(name, use_lock_file=False)
    con1._allow_reading_old_format = True
    newname = name[:-2] + 'new.db'
    with connect(newname, create_indices=False, use_lock_file=False) as con2:
        row = None
        for row in con1.select():
            kvp = row.get('key_value_pairs', {})
            con2.write(row, data=row.get('data'), **kvp)
        
        assert row is not None, 'Your database is empty!'
        
    c = con2._connect()
    for statement in index_statements:
        c.execute(statement)
    c.commit()

    os.rename(name, name[:-2] + 'old.db')
    os.rename(newname, name)
    
    
def main():
    parser = optparse.OptionParser()
    opts, args = parser.parse_args()
    for name in args:
        convert(name)
  
        
if __name__ == '__main__':
    main()
