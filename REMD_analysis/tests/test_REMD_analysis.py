"""
Unit and regression test for the REMD_analysis package.
"""

# Import package, test suite, and other packages as needed
import pytest
import sys
from argparse import Namespace
import REMD_analysis as RA

RA_test = RA.REMDAnalysis(['PLCpep7.log'])

def test_REMD_analysis_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "REMD_analysis" in sys.modules

def test_initialize():
    parser = RA.initialize(['-l', 'PLCpep7.log', '-p', 'PLCpep7'])
    assert parser == Namespace(log=['PLCpep7.log'], prefix='PLCpep7')
    assert parser.log == ['PLCpep7.log']
    assert parser.prefix == 'PLCpep7'





    

    

