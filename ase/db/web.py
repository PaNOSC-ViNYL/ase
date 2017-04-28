import os

from ase.utils import Lock


functions = []


def creates(*filenames):
    def decorator(func):
        def newfunc(row, prefix='', tmpdir='.'):
            with Lock('ase.db.web.lock'):
                for filename in filenames:
                    try:
                        os.remove(filename)
                    except FileNotFoundError:
                        pass
                func(row)
                for filename in filenames:
                    path = os.path.join(tmpdir, prefix + filename)
                    if os.path.isfile(filename):
                        os.rename(filename, path)
                    else:
                        with open(path, 'w'):
                            pass

        functions.append(newfunc)

        return newfunc

    return decorator
