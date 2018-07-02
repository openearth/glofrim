
import re, codecs, sys
from configparser import ConfigParser
from collections import OrderedDict

"""
config functios
"""
def configread(config_fn, encoding='utf-8', cf=ConfigParser()):
    """read model configuration from file"""
    cf.optionxform=str # preserve capital letter
    with codecs.open(config_fn, 'r', encoding=encoding) as fp:
        cf.read_file(fp)
    return cf 

def configwrite(config, config_fn, encoding='utf-8', cf=ConfigParser(), **kwargs):
    """write model configuration to file"""
    with codecs.open(config_fn, 'w', encoding=encoding) as fp:
        cf.write(fp)

def config2dict(cf):
    """parse config to dict"""
    out_dict = OrderedDict((sec, OrderedDict((opt, cf.get(sec, opt))
                            for opt in cf.options(sec)))
                            for sec in cf.sections())
    return out_dict