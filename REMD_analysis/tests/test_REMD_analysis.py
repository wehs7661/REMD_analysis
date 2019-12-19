"""
Unit and regression test for the REMD_analysis package.
"""

# Import package, test suite, and other packages as needed
import REMD_analysis
import pytest
import sys

def test_REMD_analysis_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "REMD_analysis" in sys.modules
