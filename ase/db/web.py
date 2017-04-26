import os

from ase.utils import Lock


functions = []


def creates(*filenames):
    def decorator(func):

        def newfunc(prefix, row):
            with Lock('ase.db.web.lock'):
                for filename in filenames:
                    os.remove(filename)
                func(row)
                for filename in filenames:
                    if os.path.isfile(filename):
                        os.rename(filename, prefix + '-' + filename)
                    else:
                        with open(prefix + '-' + filename):
                            pass

        functions.append(newfunc)

        return newfunc

    return decorator
