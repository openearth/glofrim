import os
from os.path import join
from datetime import datetime

# utils

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
