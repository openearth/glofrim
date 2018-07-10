
import re, codecs, sys
from configparser import ConfigParser
from collections import OrderedDict
import os
from subprocess import check_output, STDOUT, CalledProcessError

"""
config functions
"""
def configread(config_fn, encoding='utf-8', cf=ConfigParser()):
    """read model configuration from file"""
    cf.optionxform=str # preserve capital letter
    with codecs.open(config_fn, 'r', encoding=encoding) as fp:
        cf.read_file(fp)
    return cf 

def configwrite(cf, config_fn, encoding='utf-8', **kwargs):
    """write model configuration to file"""
    with codecs.open(config_fn, 'w', encoding=encoding) as fp:
        cf.write(fp)

def configget(config, attribute_name):
    """get value from config attribute_name (name:option)"""
    attrpath = attribute_name.split(":")
    if len(attrpath) == 2:
        return config.get(attrpath[0], attrpath[1])
    else:
        msg = "Attributes should follow the name:option convention"
        raise ValueError(msg)

def configset(config, attribute_name, attribute_value):
    """set attribute_value in config attribute_name (name:option)"""
    attrpath = attribute_name.split(":")
    if len(attrpath) == 2:
        return config.set(attrpath[0], attrpath[1], attribute_value)
    else:
        msg = "Attributes should follow the name:option convention"
        raise ValueError(msg)

def configattr(config):
    """return name:option list from config"""
    attr = []
    for sect in config.sections():
        for opt in config.options(sect):
            attr.append(sect + ":" + opt)
    return attr 

def config2dict(cf):
    """parse config to dict"""
    out_dict = OrderedDict((sec, OrderedDict((opt, cf.get(sec, opt))
                            for opt in cf.options(sec)))
                            for sec in cf.sections())
    return out_dict

def configcheck(ob, logger):
    """check if object has config attribute"""
    if not hasattr(ob, '_config'):
        msg = "run initialize_model first before using get_attribute_value"
        logger.warn(msg)
        raise Warning(msg)
    pass

"""
paths stuff
"""
def getabspath(path, root):
    if not os.path.isabs(path):
        path = os.path.normpath(os.path.join(root, path))
    return path

"""
subcall
"""
def subcall(msg, cwd='./'):
    output = check_output(msg, stderr=STDOUT, shell=True, cwd=cwd)