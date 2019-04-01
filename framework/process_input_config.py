# coding: utf-8
import json


def load_json_parameters_into_dictionary(config):
    params = read_json_files(config)
    params = params_refactoring(params)
    return params


def read_json_files(config):
    """reads a json file, outputs a dictionary"""
    with open(config) as jsn_std:
        jparams = json.load(jsn_std)

    # Remove empty strings and convert unicode characters to strings
    params = {}
    for key, val in jparams.iteritems():
        # Make sure all keys are strings
        _key = str(key)

        # ignore empty strings and comments
        if val == "" or _key == "#":
            pass
        # convert unicode values to strings
        elif isinstance(val, unicode):
            params[_key] = str(val)
        else:
            params[_key] = val

    return params


def params_refactoring(_params):
    """hack the parameter dictionary...be careful"""
    _params['wavelength'] = 1e-9 * 299792458 / _params['ms_nu']

    return _params


def setup_keyword_dictionary(prefix, dictionary):
    f = lambda x: filter(lambda a: a[0].startswith(x), dictionary.items())
    return dict([(key.split(prefix)[-1], val) for (key, val) in f(prefix)])
