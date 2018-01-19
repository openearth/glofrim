import os
from os.path import join
import re
from datetime import datetime

# utils
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
    # compile pattern: 1st group is "<key>" 2nd group is "= <old_value>"
    # matches cases with and without comments recognized "#"
    pattern = re.compile(r'(\w+)\s+(=.+)\s#?')

    # open files
    with open(fn_ini_template, 'r') as src, open(fn_ini_out, 'w') as dst:
        # loop through lines
        for line in src.readlines():
            # if the line does not match the pattern or any key, it is not changed
            line_out = line
            # replace default from template with kwarg value if found
            if pattern.match(line) is not None:
                kw, old_val = pattern.match(line).groups()[:2]
                if kw in kwargs:
                    if not isinstance(kwargs[kw], str):
                        msg = "kwargs should only contain values of type string"
                        raise ValueError(msg)
                    new_val = r'= {:s}\t'.format(kwargs[kw])
                    line_out = re.sub(old_val, new_val, line)
            # write line
            dst.write(line_out)

def determineSteps(d1, d2):
    """
    Computes numer of update steps based on start and endtime defined in PCR ini-file
    """
    d1 = datetime.strptime(d1, "%Y-%m-%d")
    d2 = datetime.strptime(d2, "%Y-%m-%d")
    return abs((d2 - d1).days)
