import functools
import numpy as np

inplace_methods = ['__iadd__', '__iand__', '__ifloordiv__', '__ilshift__',
                   '__imatmul__', '__imod__', '__imul__',
                   '__ior__', '__ipow__', '__irshift__', '__isub__',
                   '__itruediv__', '__ixor__']

forward_methods = ['__abs__', '__add__', '__and__',
                   '__bool__', '__complex__',
                   '__contains__',
                   '__divmod__', '__eq__',
                   '__float__', '__floordiv__', '__ge__',
                   '__gt__', '__hash__',
                   '__index__', '__int__',
                   '__invert__',
                   '__iter__', '__le__', '__len__',
                   '__lshift__', '__lt__', '__matmul__', '__mod__', '__mul__',
                   '__ne__',
                   '__neg__', '__or__', '__pos__', '__pow__',
                   '__radd__',
                   '__rand__', '__rdivmod__', '__reduce__', '__reduce_ex__',
                   '__rfloordiv__', '__rlshift__', '__rmatmul__',
                   '__rmod__', '__rmul__',
                   '__ror__', '__rpow__', '__rrshift__', '__rshift__',
                   '__rsub__',
                   '__rtruediv__', '__rxor__',
                   '__sub__',
                   '__truediv__', '__xor__',
                   # py2:
                   '__nonzero__']


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
    for name in forward_methods:
        if hasattr(np.ndarray, name) and not hasattr(cls, name):
            meth = forward_call(name)
            setattr(cls, name, meth)
    return cls
