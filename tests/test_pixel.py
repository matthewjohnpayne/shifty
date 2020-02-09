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
from pixel import Pixel

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

def test_WCS_methods():
    ''' Test the Pixel object and associated WCS methods '''
    print('\nWorking on test_WCS_methods() ...')
    
    
    # test basic creation of Pixel
    # -----------------------------------------

    # Need to create HDUs to populate the ImageDataSet
    T       = loader.TESSImageLoader()
    HDUs    = T._load_test_images()
    
    # Need to extract data from HDUs to form inputs to Pixel further down ...
    HDU_wcs_headers, HDU_data, HDU_unc, HDU_midtimes = T._parse_HDUs_for_ImageDataSet(HDUs)

    # Test creation of Pixel object
    P = Pixel()
    assert isinstance(P , Pixel ), \
        'IDS did not get created as expected'



    # test WCS method(s)
    # -----------------------------------------


    # test that the wcs conversion 'round-trips'
    start_pixels    = [ [1, 10, 100, 1000], [1, 10, 100, 1000] ]
    ra,dec          = np.array(WCS(HDU_wcs_headers[0]).all_pix2world(start_pixels[0], start_pixels[1], 1))
    end_pixels      = np.array(WCS(HDU_wcs_headers[0]).all_world2pix(ra, dec, 1))
    assert np.allclose(start_pixels, end_pixels), \
        'end-pixels [%r]not close enough to start pixels [%]' % (end_pixels, start_pixels)



    # test that the wcs conversion gets same result as nominal recorded in header
    '''
        part of header will look qualitatively like ...
        ...
        CRVAL1  =  37.6041907846446150 / RA at CRPIX1, CRPIX2
        CRVAL2  = -10.8342085542547150 / DEC at CRPIX1, CRPIX2
        CRPIX1  =               1045.0 / X reference pixel
        CRPIX2  =               1001.0 / Y reference pixel
        ...
    '''
    header = HDU_wcs_headers[0]
    start_pixels = [ header['CRPIX1'], header['CRPIX2'] ]
    ra,dec          = np.array(WCS(header).all_pix2world(start_pixels[0], start_pixels[1], 1))
    expectedRADEC = [ header['CRVAL1'], header['CRVAL2'] ]
    assert np.allclose( [ra,dec], expectedRADEC), \
        'end-pixels [%r]not close enough to start pixels [%r]' % (end_pixels, start_pixels)




    # use wcs to get RADEC of each pixel
    sky_coord = P.get_per_pixel_RADEC(HDU_wcs_headers[0], HDU_data[0])
    ra  = sky_coord.ra.degree
    dec = sky_coord.dec.degree
    assert shape.ra = shape.dec == (header['NAXIS2'], header['NAXIS1']), \
        'retruned arrays have an unexpected shape: shape.ra = %r ,  shape.dec = %r' % (shape.ra ,shape.dec)
    # sector4 was in the south and these test images are of camera=1, chip=1 which should be near the ecliptic
    assert np.max(dec) < 0, \
        'expected this sector4 image to be in the south ... np.max(dec)=%r' % np.max(dec)




    print(' \t Passed tests currently implemented in *test_WCS_methods()* ')







# Won't need these calls if use pytest/similar
test_WCS_methods()
