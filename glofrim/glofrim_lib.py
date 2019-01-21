
import re
import codecs
import sys
from configparser import ConfigParser
from collections import OrderedDict
import os
from os.path import dirname, basename, isfile, join
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

def configcheck(bmi, logger):
    """check if object has config attribute"""
    if not hasattr(bmi, '_config'):
        msg = "run initialize_model first before using get_attribute_value"
        logger.warn(msg)
        raise Warning(msg)
    pass

def write_config(bmi, config, config_fn, logger):
    """write adapted config to file. just before initializing
    only for models which do not allow for direct access to model config via bmi"""
    configcheck(bmi, logger)
    out_dir = dirname(config_fn)
    bname = basename(config_fn).split('.')
    new_config_fn = '{}_glofrim.{}'.format('.'.join(bname[:-1]), bname[-1])
    new_config_fn = join(out_dir, new_config_fn)
    if isfile(new_config_fn):
        os.unlink(new_config_fn)
        logger.warn("{:s} file overwritten".format(new_config_fn))
    configwrite(config, new_config_fn, encoding='utf-8')
    logger.info('Ini file written to {:s}'.format(new_config_fn))
    return new_config_fn
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

## dt check
def check_dts_divmod(dt, dt_mod):
    dt_diff = dt.total_seconds() % dt_mod.total_seconds()
    return dt_diff == 0
