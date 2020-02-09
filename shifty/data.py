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
        (i) [optional] stacked_fits_filepath
         - filepath to fits 'stack-file'
        
        
        Methods:
        --------
        observatory_position (JPL code)
        get_observatory_barycentric_positions()
        get_theta_wcs()  expose the transformation of pixel to theta space

    '''

    def __init__(self, *args) :
        
        # If an argument is supplied, try to use it to load fits 'stack-file'
        if len(args) and os.path.isfile(args[0]):
            self.load(args[0])
        
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, exc_traceback):
        # explicitly close the HDUlist object passed to __init__
        # ...
        
        # error handling ...
        if exc_type:
            print('exc_type: %r' % exc_type)
            print('exc_value: %r' % exc_value)
            print('exc_traceback: %r' % exc_traceback)
        
        # return True to handle exception...
        return True

    # -------------------------------------------------------------------------------------
    # Public Methods
    # -------------------------------------------------------------------------------------
    
    def load(self, stacked_fits_filepath ):
        '''
            Create an ImageDataSet object from a 'stack-fits-file' that was previously created.
            Stack file should contain clean data that is good-to-go
            
            inputs:
            -------
            stacked_fits_filepath
            
            returns:
            --------
            ImageDataSet object
            '''
        
        # open the fits file
        hdul = fits.open(stacked_fits_filepath)
        
        # get midtimes
        midtimes = np.array( [self._get_midtime(item.header)  for item in hdul  \
                              if isinstance(item, astropy.io.fits.hdu.image.ImageHDU )] )
    
        return hudl, midtimes


    def get_theta_coordinates(self, reference_sky_coord, HDUlist ):
        '''
            See *generate_theta_coordinates()* for a discussion of theta-coordinates
            
            Will either
            (i) read the coords from file (if previously generated)
            (ii) generate the coords & save-to-file
            
            Inputs:
            -------
            reference_sky_coord
             -
            HDUlist
             -
            
            Returns:
            --------
            theta_coords [[ ** IN WHAT FORMAT?? ** ]]
        '''
    
        # attempt to read theta-coords from HDUlist
        theta_coords = self._attempt_to_read_theta_coords_from_HDUlist(HDUlist)
        
        if theta_coords == False :
        
            # if required, do the detailed calculation of the theta coordinates
            for HDU in HDUlist :
                
                # assuming that the only ImageHDUs are the actual data
                if isinstance(HDU, astropy.io.fits.hdu.image.ImageHDU ):
                    
                    # calculate theta coordinates
                    theta_coords = self._generate_theta_coordinates( reference_sky_coord , HDU.header, HDU.data )
            
                    # append theta coordinates to HDUlist
                    self._append_theta_coordinates_to_HUDlist(theta_coords , HDUlist)
        
            # Save updated HDUlist to fits file (i.e. file now includes theta-coords)
            
        
            # Re-attempt to read theta-coords from HDUlist
            theta_coords = self._attempt_to_read_theta_coords_from_HDUlist(HDUlist)
            assert theta_coords != False, \
                'theta_coords could not be read despite attempt to generate'
    
        return theta_coords
                
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
             
            Returns:
            --------
            theta-coords
             - output of *_rot_vec_to_projection_coords()* function
             - astropy.coordinates.representation.CartesianRepresentation

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
    
    def _append_theta_coordinates_to_HUDlist( theta_coords , HDUlist ):
        '''
            Append the results of running *generate_theta_coordinates()* to HDUlist
            Includes details of the reference_sky_coord used to calculate the theta-coords
            
            Inputs:
            -------
            theta-coords
            - output of *_rot_vec_to_projection_coords()* and/or *generate_theta_coordinates()* functions
            - astropy.coordinates.representation.CartesianRepresentation
            
            HDUlist
             - HDUlist to which theta-coord data will be passed
            
            Returns:
            --------
            ???
        '''
        # Record the reference_sky_coord used to calculate the theta-coords
        # Record the original image-fits file to which these theta-coords correspond
        # Convert the theta_coords into a binary fits table
        # https://docs.astropy.org/en/stable/io/fits/#creating-a-new-fits-file
        '''
        a1 = np.array(['NGC1001', 'NGC1002', 'NGC1003'])
        a2 = np.array([11.1, 12.3, 15.2])
        hdu_table = fits.BinTableHDU.from_columns(
              [fits.Column(name='target', format='20A', array=a1),
               fits.Column(name='V_mag', format='E', array=a2)],
              name = 'Theta'
        '''
        # Append to the HDUlist
        HDUlist.append(hdu_table)
        return HDUlist

    def generate_observatory_barycentric_positions(self):
        '''
        generate barycentric-positions, one for each image
        
        we should use "wis.py"
        - I wrote it at/around TESS Ninja 2: it would be good to use it
        - N.B. wis.py currently only for space-based: would need mpc_lib functionality for ground-based
        '''
            
        # wis.py currently defaults to using heliocenter, so need to use barycenter
        return wis.Satellite(self.obscode, self.mid_times, center='BARY').posns
        

    # -------------------------------------------------------------------------------------
    # Internal Methods
    # -------------------------------------------------------------------------------------
    
    def _ensure_consistent_reference_vector(self,  *args, **kwargs ):
        '''
            Previous work has often been held-up by confusion over reference frame 
            Here we ensure that the reference vector is in a SkyCoord object
            This ensures that it can be consistently converted/interpreted
            
            inputs:
            ------
             - various possibilities allowed
             
            returns:
            --------
            SkyCoord object to be interpreted as a unit vector
            - class 'astropy.coordinates.sky_coordinate.SkyCoord'
            
        '''
    
        # for now, demand is an astropy SkyCoord object
        # - [[LATER-ON THIS CAN BE BROADENED OUT TO BE MORE FLEXIBLE IF DESIRED]]
        ref_vec = None
        for arg in args:
            if isinstance(arg, astropy.coordinates.sky_coordinate.SkyCoord):
                return arg #<<-- This doesn't seem like it would work: why indented like this ?????
    
 
    


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
            print('Problem might be because you have to pass an astropy.coordinates.sky_coordinate.SkyCoord object')
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
            print('Problem might be because you have to pass an astropy.coordinates.sky_coordinate.SkyCoord object')
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

        try:
            # convert sky-coord to cart-rep in barycentric ecliptic coords
            bary_ec_cart = reference_sky_coord.barycentricmeanecliptic.represent_as(UnitSphericalRepresentation).represent_as(CartesianRepresentation).xyz.T[0]
            # transform (rotate) using matrix from _calculate_rotation_matrix
            return vec_rep.transform( self._calculate_rotation_matrix( *bary_ec_cart ) )
        
        except Exception as error:
            print('error in _rot_vec_to_projection_coords')
            print('may be due to passing in arrays rather than (e.g.) CartesianRepresentation')
            sys.exit('%r'%error)
    




