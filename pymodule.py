import PsiMod
import re
import os
import input
import math
import warnings
from driver import *
from wrappers import *
from molutil import *
from text import *
from procutil import *


def run_fcidump(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    fcidump can be called via :py:func:`~driver.energy`.

    >>> energy('fcidump')

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    # Your plugin's PsiMod run sequence goes here
    scf_helper(name, **kwargs)
    returnvalue = PsiMod.plugin('fcidump.so')

    return returnvalue


# Integration with driver routines
procedures['energy']['fcidump'] = run_fcidump
