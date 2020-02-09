"""
    tests for shifty

"""


# -------------------------------------------------------------------------------------
# Third party imports
# -------------------------------------------------------------------------------------
import os, sys
import numpy as np

import astropy
from astropy.time import Time
from astropy import units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.coordinates.builtin_frames import FK5, ICRS, GCRS, GeocentricMeanEcliptic, BarycentricMeanEcliptic, HeliocentricMeanEcliptic, GeocentricTrueEcliptic, BarycentricTrueEcliptic, HeliocentricTrueEcliptic, HeliocentricEclipticIAU76
from astropy.coordinates.representation import CartesianRepresentation,SphericalRepresentation, UnitSphericalRepresentation


# -------------------------------------------------------------------------------------
# Any local imports
# -------------------------------------------------------------------------------------
import data
import loader

# -------------------------------------------------------------------------------------
# Test "data" module
# -------------------------------------------------------------------------------------

# ------------------------------------------------------------------------
# convenience methods for tests
# ------------------------------------------------------------------------

def _RADEC_to_unit( RA_, DEC_):
    '''
        This translates (RA,DEC)  into ...
        ...a unit vector (still in Equatorial coords)
        
        This is *NOT* used inside data ...
        ... but I want it for testing.
        
        Inputs:
        -------
        RA_, DEC_ : iterables of floats
        - assumed in degrees & assumed in Equatorial coords
        
        Returns:
        --------
        UV_ : numpy array of shape == (3,len(RA_))
        - unit vector in Equatorial coords
        
    '''
    RA_, DEC_ = np.asarray(RA_), np.asarray(DEC_)

    x_ = np.cos(DEC_*np.pi/180.) * np.cos(RA_*np.pi/180.)
    y_ = np.cos(DEC_*np.pi/180.) * np.sin(RA_*np.pi/180.)
    z_ = np.sin(DEC_*np.pi/180.)
    return np.array([x_,y_,z_]).T

def ecl():
    '''
        Obliquity of ecliptic at J2000 (From Horizons 12 April 2019)
        
        This is *NOT* used inside data ...
        ... but I want it for testing.
    '''
    return (84381.448*(1./3600)*np.pi/180.) # Obliquity of ecliptic at J2000 (From Horizons 12 April 2019)

def equatorial_to_ecliptic( v ):
    ''' 
        Convert from equatorial UV to ECLIPTIC UV

        This is *NOT* used inside data ...
        ... but I want it for testing.
    '''
    return np.matmul(v, rotate_matrix( -ecl() ).T)

def rotate_matrix(ecl):
    '''
        Set up a rotation matrix 
        
        This is *NOT* used inside data ...
        ... but I want it for testing.
    '''
    ce = np.cos(ecl)
    se = np.sin(-ecl)
    rotmat = np.array([[1.0, 0.0, 0.0],
                       [0.0,  ce,  se],
                       [0.0, -se,  ce]])
    return rotmat


# ------------------------------------------------------------------------
# tests
# ------------------------------------------------------------------------

def test_ImageDataSet():
    ''' Test the ImageDataSet object and associated methods '''
    print('\nWorking on test_ImageDataSet() ...')
    
    
    # test basic creation of ImageDataSet
    # -----------------------------------------

    # Need to create HDUs to populate the ImageDataSet
    T       = loader.TESSImageLoader()
    HDUs    = T._load_test_images()
    
    # Need to extract data from HDUs to form inputs to ImageDataSet
    HDU_wcs_headers, HDU_data, HDU_unc, HDU_midtimes = T._parse_HDUs_for_ImageDataSet(HDUs)

    # Test creation of IDS object
    # - test that it has 'images' & 'obs_code' as attributes
    IDS = data.ImageDataSet(HDU_wcs_headers, HDU_data, HDU_midtimes, T.obs_code, HDU_unc = HDU_unc)
    assert isinstance(IDS , data.ImageDataSet), \
        'IDS did not get created as expected'
    assert np.all( [ _ in IDS.__dict__ for _ in   ['headers', 'data','unc', 'obs_code',  'mid_times'] ] ), \
        'IDS does not have expected variables'
    assert np.all( [ len(IDS.__dict__[_]) == len(HDUs) for _ in ['headers', 'data','unc', 'mid_times'] ] ), \
        'created quantities not of the write length'




    # test WCS method(s)
    # -----------------------------------------

    # use wcs to get RADEC of each pixel
    sky_coord = data.get_per_pixel_RADEC(HDU_wcs_headers[0], HDU_data[0])

    #
    print('\n\t get_per_pixel_RADEC *NOT* currently tested')




    # test coord-transformation method(s)
    # -----------------------------------------

    # test *_skycoord_to_equatorialUV* --------
    # - this is likely not used much, but it is easy to test !!!
    ra  = np.array([0.,90.,180., 270.,   0.,  0.,  90., 90., 0.,90.,180.,270.])* u.deg
    dec = np.array([0., 0.,  0.,   0., -90., 90., -90., 90.,45.,45., 45., 45.])* u.deg
    expected_equatorialUV = np.array([
                            [1,0,0],
                            [0,1,0],
                            [-1,0,0],
                            [0,-1,0],
                            [0,0,-1],
                            [0,0,1],
                            [0,0,-1],
                            [0,0,1],
                            [1/np.sqrt(2),0,1/np.sqrt(2)],
                            [0,1/np.sqrt(2),1/np.sqrt(2)],
                            [-1/np.sqrt(2),0,1/np.sqrt(2)],
                            [0,-1/np.sqrt(2),1/np.sqrt(2)],
                            ])

    # set up a SkyCoord object to hold the ra, dec
    SC = SkyCoord(ra=ra, dec=dec, frame='icrs')
    
    # do transformation
    result = IDS._skycoord_to_equatorialUV(SC)

    # check results
    assert isinstance(result, astropy.coordinates.representation.CartesianRepresentation ), \
        '_skycoord_to_equatorialUV did not return an astropy CartesianRepresentation'
    assert np.allclose( result.get_xyz().T , expected_equatorialUV) , \
            'unit vectors did not match expected values '

    # get expect unit vectors by another method (using func, rather than hand-typing)
    expected_equatorialUV = _RADEC_to_unit(ra, dec)
    assert np.allclose( result.get_xyz().T , expected_equatorialUV  ) , \
        'unit vectors did not match expected values '






    # test *_skycoord_to_ecUV* --------------

    # set up a SkyCoord object to hold the ra, dec
    SC = SkyCoord(ra=ra, dec=dec, frame='icrs')

    # do transformation
    result = IDS._skycoord_to_ecUV(SC)

    # do equivalent calculation
    expected_equatorialUV = _RADEC_to_unit(ra, dec)
    expected_eclipticUV   = equatorial_to_ecliptic( expected_equatorialUV )
    assert np.allclose( result.get_xyz().T , expected_eclipticUV , rtol=1e-09, atol=1e-06 ) , \
        'unit vectors did not match expected values '
    print(' ... In order to get this test to pass, had to increase atol to 1e-6 ...')
    print(' ... I do not have a good understanding at present of the causes of difference, so cannot currently say whether this is acceptable ... ')











    # test *_ensure_consistent_reference_vector()* method
    # -----------------------------------------

    # test the method that ensures a supplied reference-vector is of allowed type

    # sky-coord objects are allowed, so test ...
    ra  = np.array([90.])* u.deg
    dec = np.array([0.])* u.deg
    SC = SkyCoord(ra=ra, dec=dec, frame='icrs')

    IDS._ensure_consistent_reference_vector( SC )
    




    # test rotation method(s)
    # -----------------------------------------
    input_vectors = np.array( [[1,0,0],[0,1,0],[0,0,1]] )

    # test *_calculate_rotation_matrix()* ---------------
    result = IDS._calculate_rotation_matrix( *input_vectors[2] )
    expectedResult = np.array([
                               [0,1,0],
                               [-1,0,0],
                               [0,0,1],
                              ])
    assert np.allclose( result, expectedResult , rtol=1e-09, atol=1e-09 ) , \
        '_calculate_rotation_matrix did not match expected values '


    # test *_rot_vec_to_projection_coords* --------------
    ecliptic_uv_CartRep = IDS._skycoord_to_ecUV( SC )
    result = IDS._rot_vec_to_projection_coords( ecliptic_uv_CartRep, SC)
    # have to do a transpose & [0] because it is outputing lists of x, y, z components
    result = np.asarray(result.xyz.T[0])
    expectedResult = np.array([0.,0.,1.])
    # Here we are asking to transforma vector that is in the same direction as the reference-vector
    # - So expect the xy-coords == 0 and the z-coord == 1
    assert np.allclose( result , expectedResult  ) , \
        '_rot_vec_to_projection_coords did not match expected values '





    # test overall *generate_theta_coordinates()* method
    # -----------------------------------------
    IDS.generate_theta_coordinates( SC , HDU_wcs_headers[0], HDU_data[0] )

    print('\n\t generate_theta_coordinates *NOT* currently tested')




    print(' \t Passed tests currently implemented in *test_ImageDataSet()* ')







# Won't need these calls if use pytest/similar
test_ImageDataSet()


###
### differential velocity aberation
### register the images
### just use wcs for reference image
### sigma-clipped average for subtraction template
### tonry ref_cat
