import warnings

class Deprecate:
    def __init__(self, obj, name, newmodule, oldmodule='ase'):
        self.obj = obj
        self.name = name
        self.newmodule = newmodule
        self.oldmodule = oldmodule

    def __call__(self, *args, **kwargs):
        message = ('%s.%s is deprecated, use %s.%s instead' %
                   (self.oldmodule, self.name, self.newmodule, self.name))
        warnings.warn(message, DeprecationWarning, stacklevel=2)
        return self.obj(*args, **kwargs)

def _dep(method):
    def _method(self, *args):
        message = ('ase.%s is deprecated, use %s.%s instead' %
                   (self.name, self.newmodule, self.name))
        warnings.warn(message, DeprecationWarning, stacklevel=2)
        return method(self, *args)
    return _method
