import os
from os.path import join
from datetime import datetime
import pdb

# utils

def determineSteps(d1, d2):
    """
    Computes numer of update steps based on start and endtime defined in PCR ini-file
    """
    d1 = datetime.strptime(d1, "%Y-%m-%d")
    d2 = datetime.strptime(d2, "%Y-%m-%d")
    return abs((d2 - d1).days)
    
    
def write_ini(fn_ini_out, fn_ini_template, **kwargs):
    """Function to fills ini-like template for all keys <kw> in kwargs.
    Note that the values in kwargs should be formatted strings and the
    keys are CASE sensitive.

    Input:
    fn_ini_out: 		new updated filename, i.e. needs to be created/defined before running the function
    fn_ini_template: 	ld filename which will serve as template for new one
    **kwarge:			whole list of

    The following lines:
    <kw> = <value> # <comment>
    <kw> =  # <comment>
    result in:
    <kw> = <new_value> # <comment>
    where: <new_value> = kwargs[<kw>]
    """
    if os.path.isfile(fn_ini_out):
        os.unlink(fn_ini_out)
    pdb.set_trace()
    # compile pattern: 1st group is "<key>" 2nd group is "= <old_value>"
    # matches cases with and without comments recognized "#", !"
    # match independent of capital letters
    kwargs = {kw.lower(): kwargs[kw] for kw in kwargs}
    # open files
    with open(fn_ini_template, 'r') as src, open(fn_ini_out, 'w') as dst:
        # loop through lines
        for line in src.readlines():
            # if the line does not match the pattern or any key, it is not changed
            line_out = line
            # replace default from template with kwarg value if found
            if '=' in line:
                kw, old_val = line.split('=')
                kw = kw.strip().lower()
                if kw in kwargs:
                    # keep comments behind ! and #
                    old_val = old_val.split('!')[0].split('#')[0].strip()
                    # replace old value
                    line_out = line.replace(old_val, kwargs[kw])
            # write line
            dst.write(line_out)

