"""These functions can be helpful to have attached to
an atoms object.
Including helper methods to enable attachment of
parametrization of the structure.
"""
import types


def get_raw_score(self):
    return self.info['key_value_pairs']['raw_score']


def set_raw_score(self, score):
    self.info['key_value_pairs']['raw_score'] = score


def enable_raw_score_methods(a):
    if 'key_value_pairs' not in a.info:
        a.info['key_value_pairs'] = {}
    a.set_raw_score = types.MethodType(set_raw_score, a)
    a.get_raw_score = types.MethodType(get_raw_score, a)


def get_neighbor_list(self):
    keys = self.info.keys()
    data = self.info['data'].keys()
    if 'data' in keys and 'neighborlist' in data:
        return self.info['data']['neighborlist']
    else:
        return None


def set_neighbor_list(self, neighbor_list):
    self.info['data']['neighborlist'] = neighbor_list


def get_parametrization(self):
    keys = self.info.keys()
    data = self.info['data'].keys()
    if 'data' in keys and 'parametrization' in data:
        return self.info['data']['parametrization']
    else:
        raise ValueError('Trying to get the parametrization before it is set!')


def set_parametrization(self, parametrization):
    self.info['data']['parametrization'] = parametrization


def enable_parametrization_methods(a):
    if 'data' not in a.info:
        a.info['data'] = {}
    a.set_neighbor_list = types.MethodType(set_neighbor_list, a)
    a.get_neighbor_list = types.MethodType(get_neighbor_list, a)
    a.set_parametrization = types.MethodType(set_parametrization, a)
    a.get_parametrization = types.MethodType(get_parametrization, a)
