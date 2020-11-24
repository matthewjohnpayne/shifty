'''
    Tests of Classes / Methods for dealing with known objects.
'''


# -----------------------------------------------------------------------------
# Third party imports
# -----------------------------------------------------------------------------
import os
import sys
import numpy as np

# -----------------------------------------------------------------------------
# Any local imports
# -----------------------------------------------------------------------------
sys.path.append(os.path.dirname(os.path.dirname(
                os.path.realpath(__file__))))
from shifty import known


# -----------------------------------------------------------------------------
# Constants and Test data
# -----------------------------------------------------------------------------

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__))), 'dev_data')


# -----------------------------------------------------------------------------
# Test "known" module
# -----------------------------------------------------------------------------


def test_empty_instantiation():
    ''' Test the class instantiation with no input'''

    # Test creation of Downloader object
    K = known.Known()
    assert isinstance(K, known.Known), \
        'Object did not get created as expected'

    print('\t completed tests of test_Downloader ')


# Won't need these calls if use pytest/similar
test_empty_instantiation()


# -------------------------------------------------------------------------
# Test data & convenience functions
# -------------------------------------------------------------------------


def _interpolate_radec(times, inputJRD):
    '''
        Interpolate the RA & Dec at the input times
        input:
        times - array of times for output
        inputJRD - tuple of JD_, RA_ and Dec_ of data for interpolation
    '''

    JD_, RA_, Dec_ = inputJRD

    # Interpolate the RA & Dec at the input times
    return np.interp(times, JD_, RA_), np.interp(times, JD_, Dec_)


def _radec_from_file(obj='sedna'):  # Thus also works for 101583
    JD_, RA_, Dec_ = np.genfromtxt(os.path.join(DATA_DIR, obj + '_ephem.txt'),
                                   usecols=(0, 1, 2), unpack=True)
    return JD_, RA_, Dec_


# End of file
