#
#@BEGIN LICENSE
#
# fcidump by James Spencer, a plugin to:
#
# PSI4: an ab initio quantum chemistry software package
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#@END LICENSE
#

import psi4
import psi4.driver.p4util as p4util
from psi4.driver.procedures import proc_util


def run_fcidump(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    fcidump can be called via :py:func:`~driver.energy`.

    >>> energy('fcidump')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = psi4.driver.scf_helper(name, **kwargs)

    proc_util.check_iwl_file_from_scf_type(psi4.core.get_option('SCF', 'SCF_TYPE'), ref_wfn)
    test_wfn = psi4.core.plugin('fcidump.so', ref_wfn)

    return test_wfn


# Integration with driver routines
psi4.driver.procedures['energy']['fcidump'] = run_fcidump
