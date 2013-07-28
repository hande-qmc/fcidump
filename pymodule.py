import psi4
import re
import os
import inputparser
import math
import warnings
from driver import *
from wrappers import *
from molutil import *
import p4util
from psiexceptions import *


def run_fcidump(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    fcidump can be called via :py:func:`~driver.energy`.

    >>> energy('fcidump')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    scf_helper(name, **kwargs)
    returnvalue = psi4.plugin('fcidump.so')

    return returnvalue

# Integration with driver routines
procedures['energy']['fcidump'] = run_fcidump
