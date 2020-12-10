'''
    Tests of Classes / Methods for dealing with known objects.
'''


# -----------------------------------------------------------------------------
# Third party imports
# -----------------------------------------------------------------------------
import os
import sys
import numpy as np
import pytest
from pytest import mark

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

all_times = np.arange(2458436.5, 2458464.5, 7)
# Some known RA and Dec values, for testing.
rd_sedna_568 = np.array([[57.13679543, 57.06643479, 56.9959375, 56.92650964],
                         [7.6546818, 7.63864473, 7.62472685, 7.61317996]])
rd_sedna_000 = np.array([[57.13678729, 57.06642135, 56.99591895, 56.92648628],
                         [7.65467493, 7.63863767, 7.62471952, 7.61317227]])
rd_101583_568 = np.array([[59.88481527, 58.23262941, 56.55372248, 54.95238514],
                          [8.51904281, 8.17400092, 7.9271755, 7.79798929]])
rd_101583_000 = np.array([[59.88447843, 58.23191998, 56.55265927, 54.9510077],
                          [8.51869996, 8.1736283, 7.92677384, 7.79756144]])
aeiOoME_sedna = [484.5211488, 0.8426141, 11.93068, 144.24763,  # a, e, i, O
                 311.35275, 358.11740, 2459000.5]  # o, M, Epoch

# -----------------------------------------------------------------------------
# Test "known" module
# -----------------------------------------------------------------------------


def test_empty_instantiation():
    '''Test the class instantiation with no input'''
    # Test creation of Known object
    K = known.Known()
    assert isinstance(K, known.Known), \
        'Object did not get created as expected'

    print('\t Completed test_empty_instantiation.')


names_of_variables = ('object_name', 'obs_code', 'times')
values_for_each_test = [
   ('Sedna', '568', all_times),
   ('101583', '568', all_times),
   ('Sedna', '000', all_times),
   ('101583', '000', all_times),
   ('Sedna', '500@-95', all_times),
   ('101583', '500@-95', all_times),
 ]
@pytest.mark.parametrize(names_of_variables, values_for_each_test)
def test__get_object_RADEC_from_horizons(object_name, obs_code, times):
    '''Test the private method for getting RA/Dec for
       a known object from Horizons.'''
    K = known.Known()
    K._get_object_RADEC_from_horizons(obs_code=obs_code, times=times,
                                      object_name=object_name)
    # Test that K.RA, K.Dec exist and have expected type.
    assert isinstance(K.RA, np.ndarray)
    assert isinstance(K.Dec, np.ndarray)
    # Test that K.RA, K.Dec have expected shape.
    assert np.shape(K.RA) == np.shape(times)
    assert np.shape(K.Dec) == np.shape(times)
    # Test that K.RA, K.Dec have expected values.
    expect_RADec = _radec_interp(times, _radec_from_file(object_name, obs_code))
    assert np.all(np.isclose(K.RA, expect_RADec[0], atol=0.00000001, rtol=0))
    assert np.all(np.isclose(K.Dec, expect_RADec[1], atol=0.00000001, rtol=0))
    print('\t Completed test__get_object_RADEC_from_horizons.')


names_of_variables = ('object_name', 'obs_code', 'times')
values_for_each_test = [
   ('Sedna', '568', all_times),
   ('101583', '568', all_times),
   ('Sedna', '000', all_times),
   ('101583', '000', all_times),
   ('Sedna', '500@-95', all_times),
   ('101583', '500@-95', all_times),
 ]
@pytest.mark.parametrize(names_of_variables, values_for_each_test)
def test_get_known_RADEC_name(object_name, obs_code, times):
    '''Test the method for getting RA/Dec for a known object, using name.'''
    K = known.Known()
    K.get_known_RADEC(obs_code=obs_code, times=times, object_name=object_name)
    # Test that K.RA, K.Dec, K.times & K.obs_code exist and have expected type.
    assert isinstance(K.RA, np.ndarray)
    assert isinstance(K.Dec, np.ndarray)
    assert isinstance(K.times, type(times))
    assert isinstance(K.obs_code, str)
    # Test that K.RA, K.Dec, K.times have expected shape.
    assert np.shape(K.RA) == np.shape(times)
    assert np.shape(K.Dec) == np.shape(times)
    assert np.shape(K.times) == np.shape(times)
    # Test that K.RA, K.Dec, K.times and K.obs_code have expected values.
    expect_RADec = _radec_interp(times, _radec_from_file(object_name, obs_code))
    assert np.all(np.isclose(K.RA, expect_RADec[0], atol=0.00000001, rtol=0))
    assert np.all(np.isclose(K.Dec, expect_RADec[1], atol=0.00000001, rtol=0))
    assert np.all(K.times == times)
    assert K.obs_code == obs_code
    print('\t Completed test_get_known_RADEC_name.')


names_of_variables = ('object_name', 'obs_code', 'times')
values_for_each_test = [
   ('Sedna', '568', all_times),
   ('101583', '568', all_times),
   ('Sedna', '000', all_times),
   ('101583', '000', all_times),
   ('Sedna', '500@-95', all_times),
   ('101583', '500@-95', all_times),
 ]
@pytest.mark.parametrize(names_of_variables, values_for_each_test)
def test_instantiate_with_object_name(object_name, obs_code, times):
    '''Test the class instantiation with an object name.'''
    K = known.Known(obs_code=obs_code, times=times, object_name=object_name)
    # Test that K.RA, K.Dec, K.times & K.obs_code exist and have expected type.
    assert isinstance(K.RA, np.ndarray)
    assert isinstance(K.Dec, np.ndarray)
    assert isinstance(K.times, type(times))
    assert isinstance(K.obs_code, str)
    # Test that K.RA, K.Dec, K.times have expected shape.
    assert np.shape(K.RA) == np.shape(times)
    assert np.shape(K.Dec) == np.shape(times)
    assert np.shape(K.times) == np.shape(times)
    # Test that K.RA, K.Dec, K.times and K.obs_code have expected values.
    expect_RADec = _radec_interp(times, _radec_from_file(object_name, obs_code))
    assert np.all(np.isclose(K.RA, expect_RADec[0], atol=0.00000001, rtol=0))
    assert np.all(np.isclose(K.Dec, expect_RADec[1], atol=0.00000001, rtol=0))
    assert np.all(K.times == times)
    assert K.obs_code == obs_code
    print('\t Completed test_instantiate_with_object_name.')


names_of_variables = ('orbit', 'obs_code', 'times')
values_for_each_test = [
    pytest.param(aeiOoME_sedna, '568', all_times,
                 marks=mark.xfail(reason='Functionality not implemented.')),
 ]
@pytest.mark.parametrize(names_of_variables, values_for_each_test)
def test__get_orbit_RADEC(orbit, obs_code, times):
    '''Test the private method for getting RA/Dec from an orbit.'''
    K = known.Known()
    K._get_orbit_RADEC(obs_code=obs_code, times=times, orbit=orbit)
    # Test that K.RA, K.Dec exist and have expected type.
    assert isinstance(K.RA, np.ndarray)
    assert isinstance(K.Dec, np.ndarray)
    # Test that K.RA, K.Dec have expected shape.
    assert np.shape(K.RA) == np.shape(times)
    assert np.shape(K.Dec) == np.shape(times)
    # Test that K.RA, K.Dec have expected values.
    expect_RADec = _radec_interp(times, _radec_from_file('Sedna', obs_code))
    assert np.all(np.isclose(K.RA, expect_RADec[0], atol=0.00000001, rtol=0))
    assert np.all(np.isclose(K.Dec, expect_RADec[1], atol=0.00000001, rtol=0))
    print('\t Completed test__get_orbit_RADEC.')


names_of_variables = ('orbit', 'obs_code', 'times')
values_for_each_test = [
    pytest.param(aeiOoME_sedna, '568', all_times,
                 marks=mark.xfail(reason='Functionality not implemented.')),
 ]
@pytest.mark.parametrize(names_of_variables, values_for_each_test)
def test_get_known_RADEC_orbit(orbit, obs_code, times):
    '''Test the method for getting RA/Dec for a known object, using orbit.'''
    K = known.Known()
    K.get_known_RADEC(obs_code=obs_code, times=times, orbit=orbit)
    # Test that K.RA, K.Dec, K.times & K.obs_code exist and have expected type.
    assert isinstance(K.RA, np.ndarray)
    assert isinstance(K.Dec, np.ndarray)
    assert isinstance(K.times, type(times))
    assert isinstance(K.obs_code, str)
    # Test that K.RA, K.Dec, K.times have expected shape.
    assert np.shape(K.RA) == np.shape(times)
    assert np.shape(K.Dec) == np.shape(times)
    assert np.shape(K.times) == np.shape(times)
    # Test that K.RA, K.Dec, K.times and K.obs_code have expected values.
    expect_RADec = _radec_interp(times, _radec_from_file('Sedna', obs_code))
    assert np.all(np.isclose(K.RA, expect_RADec[0], atol=0.00000001, rtol=0))
    assert np.all(np.isclose(K.Dec, expect_RADec[1], atol=0.00000001, rtol=0))
    assert np.all(K.times == times)
    assert K.obs_code == obs_code
    print('\t Completed test_get_known_RADEC_orbit.')


names_of_variables = ('object_name', 'obs_code', 'times', 'expected_XYZ')
values_for_each_test = [
    pytest.param('Sedna', '568', all_times, [0, 0, 0],  # xyz_sedna_568,
                 marks=mark.xfail(reason='Test not fully implemented.')),
 ]
@pytest.mark.parametrize(names_of_variables, values_for_each_test)
def test__get_object_XYZ_from_horizons(object_name, obs_code, times,
                                       expected_XYZ):
    '''Test the private method for getting RA/Dec for
       a known object from Horizons.'''
    K = known.Known()
    K._get_object_XYZ_from_horizons(obs_code=obs_code, times=times,
                                    object_name=object_name)
    # Test that K.XYZ exist and have expected type.
    assert isinstance(K.XYZ, np.ndarray)
    # Test that K.XYZ have expected shape.
    assert np.shape(K.XYZ) == np.shape([times, 3])
    # Test that K.XYZ have expected values.
    assert np.all(np.isclose(K.XYZ, expected_XYZ, atol=0.00000001, rtol=0))
    print('\t Completed test__get_object_XYZ_from_horizons.')


# -------------------------------------------------------------------------
# Test data & convenience functions
# -------------------------------------------------------------------------


def _radec_interp(times, inputJRD):
    '''
        Interpolate the RA & Dec at the input times
        input:
        times - array of times for output
        inputJRD - tuple of JD_, RA_ and Dec_ of data for interpolation
    '''

    JD_, RA_, Dec_ = inputJRD

    # Interpolate the RA & Dec at the input times
    return np.interp(times, JD_, RA_), np.interp(times, JD_, Dec_)


def _radec_from_file(obj='sedna', obs_code='C57'):  # Thus also works for 101583
    if obs_code == '500@-95':
        obs_code = 'C57'
    filename = obj + '_ephem_' + obs_code + '.txt'
    JD_, RA_, Dec_ = np.genfromtxt(os.path.join(DATA_DIR, filename),
                                   delimiter=(17, 5, 13, 13),
                                   usecols=(0, 2, 3), unpack=True)
    return JD_, RA_, Dec_


# Won't need these calls if use pytest/similar
if __name__ == '__main__':
    test_empty_instantiation()
    test_get_known_RADEC_name('Sedna', '568', all_times)
    test_instantiate_with_object_name('Sedna', '568', all_times)

# End of file
