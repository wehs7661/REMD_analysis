"""
REMD_analysis
`REMD_analysis` is a Python package of data analysis tools for replica exchange molecular dynamics (REMD) simulations.
"""

# Add imports here
from .REMD import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
