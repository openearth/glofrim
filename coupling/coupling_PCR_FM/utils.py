import os
from os.path import join
from datetime import datetime
import re, codecs
from configparser import ConfigParser
from collections import OrderedDict
from subprocess import check_output, STDOUT, CalledProcessError


# utils
def set_values_in_array(vals, idx, update_vals):
    """
    Sets new values in an existing array at specified locations. Locations are provided by a list of (y, x) indices
    :param vals: array (numpy) of values that will be updated
    :param idx: list of (y, x) indices
    :param update_vals: single value or list of values (equal size as idx) for updating
    :return: vals: updated list
    """
    vals[zip(*idx)] = update_vals
    return vals

def determineSteps(d1, d2):
    """
    Computes numer of update steps based on start and endtime defined in PCR ini-file
    """
    d1 = datetime.strptime(d1, "%Y-%m-%d")
    d2 = datetime.strptime(d2, "%Y-%m-%d")
    return abs((d2 - d1).days)


def write_ini(fn_ini_out, fn_ini_template, ignore='#', **kwargs):
    """Function to fills ini-like template for all keys <kw> in kwargs.
    Note that the values in kwargs are compared case insensitive.
    Input:
    fn_ini_out: 		new updated filename, i.e. needs to be created/defined before running the function
    fn_ini_template: 	ld filename which will serve as template for new one
    ignore:             per line, the part of the string behind the ignore token are ignored (default = '#')
    **kwargs:			whole list of key-word argmunts to be replaced
    The following lines:
    <kw> = <value> # <comment>
    <kw> =  # <comment>
    result in:
    <kw> = <new_value> # <comment>
    where: <new_value> = kwargs[<kw>]
    """
    if os.path.isfile(fn_ini_out):
        os.unlink(fn_ini_out)

    # match independent of capital letters
    kwargs = {kw.lower(): kwargs[kw] for kw in kwargs}

    # open files
    with open(fn_ini_template, 'r') as src, open(fn_ini_out, 'w') as dst:
        # loop through lines
        for line in src.readlines():
            # if the line does not match the pattern or any key, it is not changed
            line_out = line
            # split line into settings and comment part
            line_split = line.split(ignore)
            if len(line_split) == 1: # no comments
                setting = line
                comment = ''
            elif len(line_split) >= 1: # with comments
                setting = line_split[0]
                comment = ignore.join(line_split[1:])
            # replace default from template with kwarg value if found
            if '=' in setting:
                # split and strip key-word and value
                kw, old_val = setting.split('=')[:2]
                old_val = old_val.strip()
                kw = kw.strip().lower()
                if kw in kwargs:
                    if old_val == '':
                        # no old value found, only whitespaces. insert new value behind '=' token
                        line_out = setting.replace('=', '= {:s} '.format(kwargs[kw]))
                    else:
                        # replace old value with new value from kwargs
                        line_out = setting.replace(old_val, kwargs[kw])
                    # add comments back to output line
                    if len(comment) > 0:
                        line_out = '{:s}{:s}{:s}'.format(line_out, ignore, comment)
            # write line
            dst.write(line_out)

def subcall(msg, cwd='./'):
    # try:
    output = check_output(msg, stderr=STDOUT, shell=True, cwd=cwd)
    # except CalledProcessError as exc:
    #     logger.error(exc.output)

def config_to_dict(config_fn, encoding='utf-8',
                   cf=ConfigParser()):
    "read config file to dictionary"
    cf.optionxform=str # preserve capital letter
    with codecs.open(config_fn, 'r', encoding=encoding) as fp:
        cf.read_file(fp)
        out_dict = OrderedDict((sec, OrderedDict((opt, cf.get(sec, opt))
                                      for opt in cf.options(sec)))
                                for sec in cf.sections())
    return out_dict

def dict_to_config(config, config_fn, encoding='utf-8',
                   cf=ConfigParser(), **kwargs):
    "read config file to dictionary"
    if not isinstance(config, dict):
        raise ValueError("config argument should be of type dictionary")
    cf.read_dict(config)
    with codecs.open(config_fn, 'w', encoding=encoding) as fp:
        cf.write(fp)

class NamConfigParser(ConfigParser):
    def __init__(self, **kwargs):
        defaults = dict(comment_prefixes=('!', '/'),
                        inline_comment_prefixes=('!'),
                        delimiters=('='))
        defaults.update(**kwargs)
        super(ConfigParser, self).__init__(**defaults)
        self.SECTCRE = re.compile(r"&(?P<header>[^]]+)")

    def write(self, fp, space_around_delimiters=False):
        """Write an .ini-format representation of the configuration state.
        If `space_around_delimiters' is True (the default), delimiters
        between keys and values are surrounded by spaces.
        """
        super(ConfigParser, self).write(fp, space_around_delimiters=space_around_delimiters)

    def _write_section(self, fp, section_name, section_items, delimiter):
        """Write a single section to the specified `fp'."""
        fp.write(u"&{}\n".format(section_name))
        for key, value in section_items:
            value = self._interpolation.before_write(self, section_name, key,
                                                     value)
            if value is not None or not self._allow_no_value:
                value = delimiter + str(value).replace('\n', '\n\t')
            else:
                value = ""
            fp.write("{}{}\n".format(key, value))
        fp.write("/\n")
