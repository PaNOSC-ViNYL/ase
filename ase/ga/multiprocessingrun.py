""" Class for handling several simultaneous jobs.
The class has not been tested yet.
"""
from multiprocessing import Pool
import time
import copy_reg
import types
from ase.io import write, read

class MultiprocessingRun(object):

    """ Class that allows for the simultaneous relaxation of
    several candidates on the same computer.
    The method is based on starting each relaxation with an
    external python script and then monitoring when the
    relaxations are done adding in the resulting structures
    to the database.

        Parameters:
    data_connection: DataConnection object.
    tmp_folder: Folder for temporary files
    n_simul: The number of simultaneous relaxations.
    calc_script: Reference to the relaxation script.
    """
    def __init__(self, data_connection, relax_function,
                 tmp_folder, n_simul=None):
        self.dc = data_connection
        self.pool = Pool(n_simul)
        self.relax_function = relax_function
        self.tmp_folder = tmp_folder
        self.results = []

    def relax(self, a):
        self.dc.mark_as_queued(a)
        fname = '{0}/cand{1}.traj'.format(self.tmp_folder,
                                          a.info['confid'])
        write(fname, a)
        self.results.append(self.pool.apply_async(self.relax_function, [fname]))
        self._cleanup()
        
    def _cleanup(self):
        for r in self.results:
            if r.ready() and r.successful():
                fname = r.get()
                a = read(fname)
                self.dc.add_relaxed_step(a)
                self.results.remove(r)
                
    def finish_all(self):
        while len(self.results) > 0:
            self._cleanup()
            time.sleep(2.)

def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    if func_name.startswith('__') and not func_name.endswith('__'):
        #deal with mangled names
        cls_name = cls.__name__.lstrip('_')
        func_name = '_%s%s' % (cls_name, func_name)
    return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
    if obj and func_name in obj.__dict__:
        cls, obj = obj, None # if func_name is classmethod
        try:
            for cls in cls.__mro__:
                try:
                    func = cls.__dict__[func_name]
                except KeyError:
                    pass
                else:
                    break
        except AttributeError:
            func = cls.__dict__[func_name]
        return func.__get__(obj, cls)
    
# def _pickle_method(m):
#     if m.im_self is None:
#         return getattr, (m.im_class, m.im_func.func_name)
#     else:
#         return getattr, (m.im_self, m.im_func.func_name)

copy_reg.pickle(types.MethodType, _pickle_method)
