'''
Tests of Classes / Methods in the transformer module.
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
from shifty import transformer

# -----------------------------------------------------------------------------
# Constants and Test data
# -----------------------------------------------------------------------------

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__))), 'dev_data')
isq2 = 1 / (2 ** 0.5)

# -----------------------------------------------------------------------------
# Test "transformer" module
# -----------------------------------------------------------------------------

def test_empty_instantiation():
    '''
    Test the class instantiation with no input
    '''
    # Test creation of Known object
    T = transformer.Transformer()
    assert isinstance(T, transformer.Transformer), \
        'Object did not get created as expected'

    print('\t Completed test_empty_instantiation.')


names_of_variables = ('reference_radec', 'expected_projection_matrix')
values_for_each_test = [
   ([0.0, 0.0], np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])),
   ([90.0, 0.0], np.array([[-1, 0, 0], [0, 0, 1], [0, 1, 0]])),
   ([180.0, 0.0], np.array([[0, -1, 0], [0, 0, 1], [-1, 0, 0]])),
   ([270.0, 0.0], np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])),
   ([0.0, 90], np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])),
   ([0.0, -90], np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])),
   ([0.0, 45.0], np.array([[0, 1, 0], [-isq2, 0, isq2], [isq2, 0, isq2]])),
   ([45.0, 0.0], np.array([[-isq2, isq2, 0], [0, 0, 1], [isq2, isq2, 0]])),
 ]
@pytest.mark.parametrize(names_of_variables, values_for_each_test)
def test__radec_to_proj_matrix(reference_radec, expected_projection_matrix):
    '''
    Test the T._radec_to_proj_matrix method that calculates the projection
    matrix from the reference coordinates.
    '''
    T = transformer.Transformer()
    assert np.all(np.isclose(T._radec_to_proj_matrix(reference_radec),
                             expected_projection_matrix)),\
                                     'Projection matrix is wrong!'

names_of_variables = ('times', 'time0', 'radec0', 'expected_projection_matrix')
values_for_each_test = [
   (np.array([2457627.9, 2457628.9, 2457637.9, 2457727.9, 2458627.9]), 
       2457627.9, [0.0, 0.0], np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])),
   (2457627.9, 2457627.9, [90.0, 0.0],
       np.array([[-1, 0, 0], [0, 0, 1], [0, 1, 0]])),
   (2457627.9, 2457628.9, [180.0, 0.0],
       np.array([[0, -1, 0], [0, 0, 1], [-1, 0, 0]])),
 ]
@pytest.mark.parametrize(names_of_variables, values_for_each_test)
def test_instantiation_with_data(times, time0, radec0,
        expected_projection_matrix):
    '''
    Test the class instantiation with input arguments.
    '''
    T = transformer.Transformer(times=times, time0=time0, radec0=radec0)
    assert isinstance(T, transformer.Transformer), \
        'Object did not get created as expected!'
    assert np.all(np.isclose(T.proj_mat, expected_projection_matrix)),\
            'Projection matrix is wrong!'
    assert np.all(T.times==times), 'times not assigned correctly!'
    assert T.time0==time0, 'reference time not recorded correctly!'

    print('\t Completed test_instantiation_with_data.')


names_of_variables = ('times', 'obs_code')
values_for_each_test = [
    (np.array([2457627.9, 2457628.9, 2457637.9, 2457727.9, 2458627.9]), '568'),
    (2457627.86493123, '568'),
    (2457627.86493123, '500'),
    (2457627.86493123, 'I11'),
 ]
@pytest.mark.parametrize(names_of_variables, values_for_each_test)
def test_get_heliocentric_equatorial_XYZ(times, obs_code):
    '''
    Test the get_heliocentic_equatorial_XYZ_from_JPL and the
    get_heliocentric_equatorial_from_MPC against each other.
    They should produce comparable results. 
    '''
    args = {'times': times, 'obs_code': obs_code}
    heMPC = transformer.get_heliocentric_equatorial_XYZ_from_MPC(**args)
    heJPL = transformer.get_heliocentic_equatorial_XYZ_from_JPL(**args)
    print(heMPC)
    print(heJPL)
    print(np.isclose(heMPC, heJPL))
    assert np.all(np.isclose(heMPC, heJPL)),\
            'JPL and MPC disagree on observer location!'

names_of_variables = ('observatory_position', 'reference_radec',
                      'expected_projected_position')
values_for_each_test = [
   ([[- 1, 1, 0]], [0, 0], [[1, 0, -1]]),
   ([[- 1, 1, 0]], [45, 0], [[2**0.5, 0, 0]]),
   ([[isq2 - 1, isq2, 0]], [0, 0], [[isq2, 0, isq2 - 1]]),
   ([[0, 2, 0]], [0, 45], [[2, 0, 0]]),
   ([[2, 0, 0]], [0, 45], [[0, -2**0.5, 2**0.5]]),
   ([[0, 0, 2]], [0, 45], [[0, 2 * isq2, 2 * isq2]]),
   ([[1, 0, 0]], [45, 45], [[-isq2, -isq2 * isq2, isq2 * isq2]]),
   ([[0, 0, 0], [-1, +1, 0], [-2, 0, 0], [-1, -1, 0]], [0, 0],
    [[0, 0, 0], [+1, 0, -1], [0, 0, -2], [-1, 0, -1]]),
 ]
@pytest.mark.parametrize(names_of_variables, values_for_each_test)
def test__do_transformation(observatory_position, reference_radec,
                            expected_projected_position):
    T = transformer.Transformer(radec0=reference_radec)
    projected_position = T._do_transformation(observatory_position)
    print(projected_position, expected_projected_position)
    print(np.isclose(projected_position, expected_projected_position))
    print(np.all(np.isclose(projected_position, expected_projected_position)))
    assert np.all(np.isclose(projected_position, expected_projected_position)),\
           "Projection/rotation wasn't correct"

#def test__xyz_observer():
#    T = transformer.Transformer(times=[2459294.187934028],
#                                time0=2459294.187934028,
#                                radec0=[0, 0])

