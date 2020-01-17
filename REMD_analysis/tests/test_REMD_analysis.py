"""
Unit and regression test for the REMD_analysis package.
"""

# Import package, test suite, and other packages as needed
import pytest
import sys
import copy
import numpy as np
from argparse import Namespace
import REMD_analysis as RA

RA_test = RA.REMDAnalysis(['test1.log'])

def test_REMD_analysis_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "REMD_analysis" in sys.modules

def test_initialize():
    parser = RA.initialize(['-l', 'PLCpep7.log', '-p', 'PLCpep7'])
    assert parser == Namespace(log=['PLCpep7.log'], prefix='PLCpep7')
    assert parser.log == ['PLCpep7.log']
    assert parser.prefix == 'PLCpep7'

def test_init():
    # this is the unit test for the __init__ function for both classes
    expected = {
        'N_states': 4,
        'dt': 0.002,
        'final_t': None,
        'finish': None,
        'n_ex': 0,
        'nex': True,
        'replex': 100.0,
        'sample_all': None,
        'start': 576
    }
    calculated = vars(RA_test)

    assert calculated == expected


def test_get_replica_data():
    time, state, t_matrix = RA_test.get_replica_data(['test1.log'])
    t_matrix_expected = np.array([[0.2485, 0.2435, 0.2605, 0.2475],
       [0.2435, 0.2766, 0.2365, 0.2435],
       [0.2605, 0.2365, 0.2465, 0.2565],
       [0.2475, 0.2435, 0.2565, 0.2525]])
    
    assert len(state) == RA_test.N_states
    assert (t_matrix == t_matrix_expected).all()

def test_plot_state_data():
    pass

def test_plot_matrix():
    pass

    
    






    

    

