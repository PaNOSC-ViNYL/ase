"""These functions can be helpful to have attached to
an atoms object.
Including helper methods to enable attachment of
parametrization of the structure.
"""
import types
import warnings

m = "enable_raw_score_methods, atoms.get_raw_score and atoms.set_raw_score "
m += "methods are deprecated and will be removed soon.\n"
m += "Use instead:\n"
m += "atoms.info['key_value_pairs']['raw_score'] directly to set and get\n"
m += "OR:\n"
m += "from ase.ga import get_raw_score, set_raw_score\n"
m += "set_raw_score(atoms, raw_score)  # to set and \n"
m += "get_raw_score(atoms)  # to get."
raw_score_message = (m)


def get_raw_score(self):
    warnings.warn(raw_score_message, FutureWarning)
    return self.info['key_value_pairs']['raw_score']


def set_raw_score(self, score):
    warnings.warn(raw_score_message, FutureWarning)
    self.info['key_value_pairs']['raw_score'] = score


def enable_raw_score_methods(a):
    warnings.warn(raw_score_message, FutureWarning)
    if 'key_value_pairs' not in a.info:
        a.info['key_value_pairs'] = {}
    a.set_raw_score = types.MethodType(set_raw_score, a)
    a.get_raw_score = types.MethodType(get_raw_score, a)


m = "enable_parametrization_methods, atoms.get_ and atoms.set_neighbor_list "
m += "and atoms.get_ and atoms.set_parametrization are deprecated and will "
m += "be removed soon, use instead:\n"
m += "atoms.info['data']['parametrization'] and "
m += "atoms.info['data']['neighborlist'] to set and get directly\n"
m += "OR:\n"
m += "from ase.ga import set_parametrization, set_neighbor_list\n"
m += "from ase.ga import get_parametrization, get_neighbor_list\n"
m += "set_parametrization(atoms, parametrization)  # to set\n"
m += "get_parametrization(atoms)  # to get"
message = (m)


def get_neighbor_list(self):
    warnings.warn(message, FutureWarning)
    keys = self.info.keys()
    data = self.info['data'].keys()
    if 'data' in keys and 'neighborlist' in data:
        return self.info['data']['neighborlist']
    else:
        return None


def set_neighbor_list(self, neighbor_list):
    warnings.warn(message, FutureWarning)
    self.info['data']['neighborlist'] = neighbor_list


def get_parametrization(self):
    warnings.warn(message, FutureWarning)
    keys = self.info.keys()
    data = self.info['data'].keys()
    if 'data' in keys and 'parametrization' in data:
        return self.info['data']['parametrization']
    else:
        raise ValueError('Trying to get the parametrization before it is set!')


def set_parametrization(self, parametrization):
    warnings.warn(message, FutureWarning)
    self.info['data']['parametrization'] = parametrization


def enable_parametrization_methods(a):
    warnings.warn(message, FutureWarning)
    if 'data' not in a.info:
        a.info['data'] = {}
    a.set_neighbor_list = types.MethodType(set_neighbor_list, a)
    a.get_neighbor_list = types.MethodType(get_neighbor_list, a)
    a.set_parametrization = types.MethodType(set_parametrization, a)
    a.get_parametrization = types.MethodType(get_parametrization, a)
