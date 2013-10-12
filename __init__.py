"""Plugin docstring.

"""
__version__ = '0.1'
__author__  = 'Psi4 Developer'

# Load Python modules
# Something weird with boost_python, but the explicit import doesn't work on
# python 2 (at least, not on CMTH machines) whereas the relative import doesn't
# work on python 3 (i.e. my laptop).
try:
    from fcidump.pymodule import *
except ImportError:
    from pymodule import *

# Load C++ plugin
import os
import psi4
plugdir = os.path.split(os.path.abspath(__file__))[0]
sofile = plugdir + '/' + os.path.split(plugdir)[1] + '.so'
psi4.plugin_load(sofile)

