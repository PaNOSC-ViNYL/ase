import sys
import numpy as np

                   #'__iand__',
                   #'__ifloordiv__',
                   #'__ilshift__',
                   #'__imod__',
                   #'__ior__',
                   #'__irshift__',
                   #'__ixor__'

py3 = sys.version_info[0] == 3

inplace_methods = ['__iadd__',
                   '__imul__',
                   '__ipow__',
                   '__isub__',
                   '__itruediv__']
if py3:
    inplace_methods.append('__imatmul__')
else:
    inplace_methods.append('__idiv__')

 #'__and__',
 #'__bool__',
 #'__complex__',
 #'__divmod__',
                   #'__float__', '__floordiv__',
                   #'__index__', '__int__',
                   #'__invert__',
                   #'__lshift__',
#'__mod__',
#'__or__',
                   #'__ror__',
 #'__rrshift__', '__rshift__',
 #'__rand__', '__rdivmod__',
 #'__reduce__', '__reduce_ex__',
 #'__rfloordiv__', '__rlshift__',
 #'__rmod__',


forward_methods = ['__abs__',
                   '__add__',
                   '__contains__',
                   '__eq__',
                   '__ge__',
                   '__getitem__',
                   '__gt__', '__hash__',
                   '__iter__', '__le__', '__len__',
                   '__lt__',
                   '__mul__',
                   '__ne__',
                   '__neg__',
                   '__pos__', '__pow__',
                   '__radd__',
                   '__rmul__',
                   '__rpow__',
                   '__rsub__',
                   '__rtruediv__', #'__rxor__',
                   '__setitem__', '__sub__',
                   '__truediv__', #'__xor__',
                   # py2:
                   #'__nonzero__'
]

if py3:
    forward_methods += ['__matmul__', '__rmatmul__']
else:
    forward_methods += ['__div__', '__rdiv__']

forward_methods += ['all', 'any', 'diagonal', 'dot', 'sum', 'ravel', 'tolist',
                    'transpose', 'tofile', 'tobytes', 'tostring']


#'argmax', 'argmin',
#'argpartition', 'argsort', 'astype',
#'dtype',
#'ndim', # ?
#'nonzero',
#'swapaxes',
#'size',

#forward_things = ['T', 'all', 'any',
#                  'diagonal', 'dot',
#                  'flat', 'flatten',
#                  'max', 'min',
#                  'ravel', 'shape',
#                  'sum',
#                  'tolist', 'transpose']


def forward_inplace_call(name):
    arraymeth = getattr(np.ndarray, name)
    def f(self, obj):
        a = self.__array__()
        arraymeth(a, obj)
        return self
    # use update_wrapper()?
    f.__name__ = name
    f.__qualname__ = name
    return f


def forward_call(name):
    arraymeth = getattr(np.ndarray, name)
    def f(self, *args, **kwargs):
        a = self.__array__()
        return arraymeth(a, *args, **kwargs)
    f.__name__ = name
    f.__qualname__ = name
    return f


def arraylike(cls):
    for name in inplace_methods:
        if hasattr(np.ndarray, name) and not hasattr(cls, name):
            meth = forward_inplace_call(name)
        setattr(cls, name, meth)
    for name in forward_methods: # + forward_things:
        assert hasattr(np.ndarray, name), name
        if hasattr(np.ndarray, name) and not hasattr(cls, name):
            meth = forward_call(name)
            setattr(cls, name, meth)
    return cls
