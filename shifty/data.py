'''
   Definition of shifty.py's fundamental data-structure
   - The "ImageDataSet"
'''

# -------------------------------------------------------------------------------------
# Third party imports
# -------------------------------------------------------------------------------------
import os, sys
from collections import OrderedDict
import numpy as np
import functools

import astropy
from astropy.time import Time
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.coordinates.matrix_utilities import rotation_matrix
from astropy.coordinates.builtin_frames import FK5, ICRS, GCRS, GeocentricMeanEcliptic, BarycentricMeanEcliptic, HeliocentricMeanEcliptic, GeocentricTrueEcliptic, BarycentricTrueEcliptic, HeliocentricTrueEcliptic, HeliocentricEclipticIAU76
from astropy.coordinates.representation import CartesianRepresentation,SphericalRepresentation, UnitSphericalRepresentation


# -------------------------------------------------------------------------------------
# Any local imports
# -------------------------------------------------------------------------------------
# N/A


# -------------------------------------------------------------------------------------
# Various class definitions for *Data Handling* in shifty
# -------------------------------------------------------------------------------------


class ImageDataSet():
    '''
        A set of images w/ times & wcs which are suitable for stacking 
        (e.g. stars have been masked, bad cadences removed, …)
        
        Inputs:
        -------
        (i) images:
        list of `astropy.io.fits.hdu.hdulist.HDUList` objects
        
        (ii) obs_code:
        string
         - mpcObs code that uniquely specifies observatory, 
         - allows determination of observatory position as a func of time
        
        Methods:
        --------
        observatory_position (JPL code)
        get_observatory_barycentric_positions()
        get_theta_wcs()  expose the transformation of pixel to theta space

    '''

    def __init__(self, HDU_wcs_headers, HDU_data, HDU_midtimes , obs_code, HDU_unc = None) :
        
        # [[ Being indicisive about how to structure ]]
        # [[ Worried about FITS headers & WCS coeffs ]]
        # [[ E.g. Believe PS1 fits/smf have 2-levels of WCS headers ]]
        # [[ So how to represent this generally ? ]]
        # [[ Will ignore the issue for now ]]

        # header part of HDUs from fits: needs to contain wcs
        self.headers = HDU_wcs_headers
        
        # data parts of HDUs from fits
        self.data       = np.asarray(HDU_data)
        
        # uncertainties parts of HDUs from fits (if any)
        self.unc        = np.asarray(HDU_unc)

        # Exposure mid-times (BJD)
        self.mid_times = HDU_midtimes
        
        # mpcObs code that uniquely specifies observatory
        self.obs_code  = obs_code

    # -------------------------------------------------------------------------------------
    # Public Methods
    # -------------------------------------------------------------------------------------

    def generate_observatory_barycentric_positions(self):
        ''' 
            generate barycentric-positions, one for each image
            
            we should use "wis.py"
             - I wrote it at/around TESS Ninja 2: it would be good to use it
             - N.B. wis.py currently only for space-based: would need mpc_lib functionality for ground-based
        '''
        
        # wis.py currently defaults to using heliocenter, so need to use barycenter
        return wis.Satellite(self.obscode, self.mid_times, center='BARY').posns

    def generate_theta_coordinates(self , reference_sky_coord , header, data ):
        ''' 
            Each pixel in an image corresponds to some RA,Dec coord (can get from wcs)
             - N.B. an (RA,Dec) can be represented as a unit-vector, u=(u_x, u_y, u_z)
             
            We want to rotate/transform the coordinates of each pixel in each image into a single, standard, set of projection-plane coordinates
        
            We are *NOT* touching flux values
            
            For each image this will return a 2D array of "theta-coordinates", one for each pixel in the image
            
            We need to think about whether/how the tangent plane is pixelized
            
            inputs:
            -------
            reference_sky_coord
            - astropy.coordinates.sky_coordinate.SkyCoord object
            - used to define a vector, in a coord frame, that will define a rotation-transformation

            header
             - 
             
            data
             -

        '''
        # do something to ensure that the reference-vector is in a consistent format &/or frame
        reference_sky_coord = self._ensure_consistent_reference_vector( reference_sky_coord )
        
        # get RA,Dec for each pixel (using wcs)
        sky_coord = get_per_pixel_RADEC(header, data)
    
        # convert ra,dec to unit-vector in ECLIPTIC coordinates
        ecliptic_uv_CartRep = self._skycoord_to_ecUV( sky_coord )
    
        # rotate ecliptic unit vector to projection plane (theta) coords
        return self._rot_vec_to_projection_coords(ecliptic_uv_CartRep ,
                                                  reference_sky_coord )

    # -------------------------------------------------------------------------------------
    # Internal Methods
    # -------------------------------------------------------------------------------------
    
    def _ensure_consistent_reference_vector(self,  *args, **kwargs ):
        '''
            Previous work has often been held-up by confusion over reference frame 
            Going to ensure that whatever is passed is made into ECLIPTIC coordinates
            
            inputs:
            ------
             - various possibilities allowed
             
            returns:
            --------
            ecliptic unit vector
            - class 'astropy.coordinates.representation.CartesianRepresentation'
            
        '''
    
        # for now, demand is an astropy SkyCoord object
        # - [[LATER-ON THIS CAN BE BROADENED OUT TO BE MORE FLEXIBLE IF DESIRED]]
        ref_vec = None
        for arg in args:
            if isinstance(arg, astropy.coordinates.sky_coordinate.SkyCoord):
                return arg
    
 
    


    def _skycoord_to_ecUV(self, sky_coord ):
        '''
            convert sky_coord (ra,dec) to unit-vector in BARYCENTRIC ECLIPTIC coordinates
            https://docs.astropy.org/en/stable/coordinates/
            https://docs.astropy.org/en/stable/coordinates/#built-in-frame-classes
            https://docs.astropy.org/en/stable/coordinates/representations.html
            https://docs.astropy.org/en/stable/api/astropy.coordinates.UnitSphericalRepresentation.html#astropy.coordinates.UnitSphericalRepresentation
         
            inputs:
            -------
            sky_coord
            - astropy.coordinates.sky_coordinate.SkyCoord object
            
            returns:
            --------
            ecliptic unit vector
             - class 'astropy.coordinates.representation.CartesianRepresentation'
        '''
        
        # Convert to barycentric, then represent as spherical-coords, then convert to unit-cartesian
        try:
            return sky_coord.barycentricmeanecliptic.represent_as(UnitSphericalRepresentation).represent_as(CartesianRepresentation)
        except Exception as error:
            print('you have to pass an astropy.coordinates.sky_coordinate.SkyCoord object')
            print(error)
            return None
    
    def _skycoord_to_equatorialUV(self, sky_coord ):
        '''
            convert  skycoord ( ra,dec) to unit-vector in ECLIPTIC coordinates
            https://docs.astropy.org/en/stable/coordinates/
            https://docs.astropy.org/en/stable/coordinates/#built-in-frame-classes
            https://docs.astropy.org/en/stable/coordinates/representations.html
            https://docs.astropy.org/en/stable/api/astropy.coordinates.UnitSphericalRepresentation.html#astropy.coordinates.UnitSphericalRepresentation
            
            [[LIKELY NOT USED IN PRACTICE : JUST PROVIDING FOR CONVENIENCE ]]
            
            inputs:
            -------
            sky_coord
            - astropy.coordinates.sky_coordinate.SkyCoord object
            
            returns:
            --------
            equatorial unit vector(s)
            - class 'astropy.coordinates.representation.CartesianRepresentation'
            
        '''
        # Convert to equatorial (icrs), then represent as spherical-coords, then convert to unit-cartesian
        try:
            return sky_coord.icrs.represent_as(UnitSphericalRepresentation).represent_as(CartesianRepresentation)
        except Exception as error:
            print('you have to pass an astropy.coordinates.sky_coordinate.SkyCoord object')
            print(error)
            return None

    @functools.lru_cache(maxsize=10)
    def _calculate_rotation_matrix(self, x_ref, y_ref, z_ref):
        ''' This routine returns the 3-D rotation matrix for the
            given reference vector.
            
            N.B. Deliberately asking for components to facilitate caching
            
            rot_mat is a rotation matrix that converts from (e,g. ecliptic)
            vectors to the projection coordinate system.
            
            The projection coordinate system has z outward,
            x parallel to increasing ecliptic longitude, and
            y northward, making a right-handed system
            
            [[ There is some possibility that this function could be replaced by astropy functionality ]]
            [[ But that would require me to understand how to use ... ]]
            [[ ... from astropy.coordinates.matrix_utilities import rotation_matrix, matrix_product ]]
            [[ rotation = rotation_matrix(45 * u.deg, axis='z') ]]
            
        '''
        
        x_ref, y_ref, z_ref = float(x_ref), float(y_ref), float(z_ref)
        r = np.sqrt(x_ref*x_ref + y_ref*y_ref + z_ref*z_ref)
        lon0 = np.arctan2(y_ref, x_ref)
        lat0 = np.arcsin(z_ref/r)
        slon0 = np.sin(lon0)
        clon0 = np.cos(lon0)
        slat0 = np.sin(lat0)
        clat0 = np.cos(lat0)
        
        rot_mat = np.array([[-slon0, clon0, 0],
                            [-clon0*slat0, -slon0*slat0, clat0],
                            [clon0*clat0, slon0*clat0, slat0 ]])
        return rot_mat
    
    
    def _rot_vec_to_projection_coords(self, vec_rep, reference_sky_coord ):
        '''
            rotate to projection plane
            https://docs.astropy.org/en/stable/api/astropy.coordinates.CartesianRepresentation.html
            https://github.com/astropy/astropy/blob/master/astropy/coordinates/matrix_utilities.py
            
            "vec_rep" and "reference_vec" need to be in consistent coordinates
            - typically we will be working in ECLIPTIC coords, but any consistent choice would work
            
            inputs:
            -------
            vec_rep
            - astropy.coordinates.representation.CartesianRepresentation
            - vectors that are to be transformed
            
            reference_sky_coord
            - astropy.coordinates.sky_coordinate.SkyCoord object
            - used to define a vector, in a coord frame, that will define a rotation-transformation
            
            returns:
            --------
            transformed_vectors
            - astropy.coordinates.representation.CartesianRepresentation
            
            
        '''

        if 1:
            # convert sky-coord to cart-rep in barycentric ecliptic coords
            bary_ec_cart = reference_sky_coord.barycentricmeanecliptic.represent_as(UnitSphericalRepresentation).represent_as(CartesianRepresentation).xyz.T[0]
            # transform (rotate) using matrix from _calculate_rotation_matrix
            return vec_rep.transform( self._calculate_rotation_matrix( *bary_ec_cart ) )
        
        else:#except Exception as error:
            print('error in _rot_vec_to_projection_coords')
            print('may be due to passing in arrays rather than (e.g.) CartesianRepresentation')
            sys.exit('%r'%error)
    




