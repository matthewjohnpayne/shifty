'''
   Contains the definition of shifty.py's fundamental data-structure
   - The "ImageDataSet"
   
   At present this contains various methods related to the
   cleaning of data.
    - MJP 2020-10-12
    - I think the cleaning-stuff should probably be moved to a completely different
      class / module
    - Perhaps make some kind of "clean" class to do the various operations
      required to from noisy fits data to cleaned "ImageDataSet"
    - I.e. by definition an "ImageDataSet" should be cleaned & read to stack
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
        Holds a set of imaging data to stack.
        
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
        
        # If arguments are supplied, try to use them to load fits 'stack-file'
        if len(args):
            self.load(args)

    def __len__(self):
        return self.data.shape[0]

    """
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
    """

    # -------------------------------------------------------------------------------------
    # Coordinate Transformation Methods
    # - *NO* WCS/reproject methods set-up as yet 
    # -------------------------------------------------------------------------------------

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
    




    # -------------------------------------------------------------------------------------
    # The methods below are for the "CLEANING" of TESS data
    # -------------------------------------------------------------------------------------



    def _clean_data(self, header , imageData, **kwargs):
        '''
            Wrapper around the various "cleaning" methods below
            - Simply provided as a means to enable simple toggling on/off of functions
            - Uses kwargs to control what is/isn't evaluated
        '''
        # dict to hold key:function mappings
        # - Have used ORDERED because at some point it might be important what order the funcs are used
        cleaning_function_dict = OrderedDict([
                              ('mask'      , self._mask_stars ),
                              ('subtract'  , self._subtract_stars ),
                              ('clip'      , self._clip_peaks ),
                              ('bad_cad'   , self._remove_bad_cadences ),
                              ('scat'      , self._remove_scattered_light_problem_areas ),
                              ('strap'     , self._remove_strap_regions ),
                              ])

        # loop over possible funcs (in order of dict)
        # [[ N.B. HDU_data should be modified in-place by functions ]]
        for key, func_to_run in cleaning_function_dict.items():
            # run a function if it is included as True in cleaning_parameters (e.g. {'mask':True})
            if key in kwargs and kwargs[key]:
                imageData = func_to_run(header , imageData, **kwargs)

        return True



    def _mask_stars(self, header , imageData, **kwargs):
        '''
            We want to remove stars in some way
            Barentsen & Payne discussed a simple mask: i.e. set pixels that contain stars to NaN
            This would be done based on GAIA positions
            
            This is *NOT* subtraction (see _subtract_stars below )
            
            Presumably only one of _mask_stars / _subtract_stars / _clip_peaks is required
            
            '''
                
                # Provide a means to only do the refcat search once
                # If NO useful refcat dictionary supplied, do search, otherwise use supplied dictionary
                # - using this implicitly assumes that all of the images are closely aligned (v. similar ra,dec ranges)
                # -
                if 'refcat_dict' not in kwargs or kwargs['refcat_dict'] == {} :
                print('in calc loop')
                kwargs['refcat_dict'] = {}
                ra,dec,pix,int_pix = RefCat().find_all_stars_on_image(header , imageData)
                kwargs['refcat_dict']['ra'], kwargs['refcat_dict']['dec'] , kwargs['refcat_dict']['pix'] , kwargs['refcat_dict']['int_pix'] = ra,dec,pix, int_pix
                    
                    
                    # Need to do something about deciding how big a mask to use, based on the source magnitude
                    # Perhaps something from photutils
                    # https://photutils.readthedocs.io/en/stable/psf.html
                    # http://docs.astropy.org/en/stable/api/astropy.convolution.discretize_model.html
                    # Perhaps using the downloaded prf
                    print(' ** WARNING: just outputting a single central mask pixel at present ** ')
                        
                        
                        # mask all of the stars
                        # - N.B. this is likely to fail if nPixels > 0 in RefCat().find_all_stars_on_image()
                        # - N.B. this alters imageData-in-place ...
                        rows, cols = kwargs['refcat_dict']['int_pix'][1] , kwargs['refcat_dict']['int_pix'][0]
                            imageData[rows, cols] = 0


        # Pre-compute the low pass filter
        fc = 0.12
        b = 0.08
        N = int(np.ceil((4 / b)))
        if not N % 2:
        N += 1
        n = np.arange(N)
        sinc_func = np.sinc(2 * fc * (n - (N - 1) / 2.0))
        window = 0.42 - 0.5 * np.cos(2 * np.pi * n / (N - 1))
        window += 0.08 * np.cos(4 * np.pi * n / (N - 1))
        sinc_func = sinc_func * window
        sinc_func = sinc_func / np.sum(sinc_func)

        def lowpass(vec, axis=-1):
        new_signal = np.convolve(vec, sinc_func)
        lns = len(new_signal)
        diff = int(
        np.abs(lns - len(vec)) / 2
        )  # low pass filter to remove high order modes
        return new_signal[diff:-diff]


        def do_pca(i, j, data, factor, q):
        xvals = np.dot(factor, data[i, j, :][q])
        return xvals


        def calc_2dbkg(flux, qual, time, fast=True):
            ''' This from eleanor '''
        flux = np.ascontiguousarray(np.moveaxis(flux, 0, -1), dtype=np.float64)
        q = qual == 0

        # build a single frame in shape of detector.
        # This was once the median image, not sure why it was changed
        med = np.percentile(flux[:, :, :], 1, axis=(2))

        # subtract off median to remove stable background emission, which is
        # especially apparent in corners
        med = med - medfilt2d(med, 25)
        print('med.shape', med.shape)

        # mask which should separate pixels dominated by starlight from background
        g = np.ma.masked_where(med < np.percentile(med, 70.0), med)
        print('g.shape', g.shape)

        modes = 21
        X = np.ascontiguousarray(flux[g.mask][:, q])
        np.random.seed(42)
        pca = PCA(n_components=modes, svd_solver="randomized")
        pca.fit(X)
        pv = pca.components_[0:modes]
        pvT = np.ascontiguousarray(pv.T)
        print('pvT.shape', pvT.shape)

        rolls = list(range(-15, 15, 6))
        N = pv.shape[1]
        vv = np.zeros((N, modes * (1 + len(rolls)) + 1), dtype=np.float64)
        vv[:, :modes] = pvT
        vv[:, -1] = 1.0
        for j, i in enumerate(rolls):
        if i != 0:
        if i > 0:
        vv[i:, (j + 1) * modes : (j + 2) * modes] = pvT[:-i]
        else:
        vv[:i, (j + 1) * modes : (j + 2) * modes] = pvT[-i:]

        vvT = vv.T
        factor = np.linalg.solve(np.dot(vvT, vv), vvT)

        nx, ny = flux.shape[:2]
        nv = vv.shape[1]
        data = np.empty((nx, ny, nv))
        data[:] = np.nan
        data[g.mask] = np.dot(X, factor.T, out=data[g.mask])

        outmeasure = np.nansum(
        np.abs(data - np.nanmean(data, axis=(0, 1)))
        / np.nanstd(data, axis=(0, 1)),
        axis=-1,
        )

        # these eigenvectors give a bad fit when something astrophysical happens
        # to these pixels, like a bright asteroid crossing the FOV. Isolate these
        # pixels and add them into the mask of pixels we don't use
        mask = ~g.mask

        metric = outmeasure / data.shape[-1]
        mask[metric > 1] = True

        x = np.arange(0, data.shape[1])
        y = np.arange(0, data.shape[0])

        xx, yy = np.meshgrid(x, y)

        count = np.zeros(data.shape[:2], dtype=np.int32)
        for i in range(np.shape(vv)[1]):
        fill_grid.fill_grid(data[:, :, i], ~mask, count)

        bkg_arr = np.zeros_like(flux)
        bkg_arr[:, :, q] = np.dot(data, vv.T)
        for i in range(nx):
        for j in range(ny):
        bkg_arr[i, j, q] = lowpass(bkg_arr[i, j, q])

        f = interp1d(
        time[q],
        bkg_arr[:, :, q],
        kind="linear",
        axis=2,
        bounds_error=False,
        fill_value="extrapolate",
        )  # then interpolate the background smoothly across
        bkg_arr[:, :, ~q] = f(
        time[~q]
        )  # all times, where we've only been learning with the good quality cadences

        return np.asfortranarray(np.moveaxis(bkg_arr, -1, 0), dtype=np.float32) , mask 


    def _subtract_stars(self, header , imageData, **kwargs):
        '''
            We want to remove stars in some way
            Holman & Payne have generally assumed some form of subtraction
            
            This is *NOT* masking (see _mask_stars )
            
            Presumably only one of _mask_stars / _subtract_stars / _clip_peaks is required
            - but I am 100% certain that Holman will at least want to experiment with subtraction
            
            Input:
            --------
            
            Returns:
            --------
            '''
                # Naive
                # - Just do subtraction of first image as template, with little/no registering (that's what Oelkers ended up doing in TESS RNAAS)
                
                # There is the DIA tool by Oelkers
                # https://iopscience.iop.org/article/10.3847/1538-3881/aad68e/meta
                # https://iopscience.iop.org/article/10.3847/1538-3881/aad68e/pdf
                # https://github.com/ryanoelkers/DIA
                
                
                
                # There is the HOTPANTS tool by Becker
                # http://web.ipac.caltech.edu/staff/fmasci/home/astro_refs/HOTPANTSsw2011.pdf
                # https://github.com/acbecker/hotpants
                
                
                pass
            
            def _clip_peaks(self, header , imageData, **kwargs):
            '''
                We want to remove stars in some way
                Could use some form of peak-clipping to remove any/all points above some ~background level
                
                Presumably only one of _mask_stars / _subtract_stars / _clip_peaks is required
                
                Input:
                --------
                
                Returns:
                --------
                '''


    def _remove_bad_cadences(self,header , imageData, **kwargs):
        '''
            In many cases it may be most straightforward to simply eliminate entire exposures
            E.g. If they have terribly high stray-light across the entire exposure
            
            Input:
            --------
            list HDUs
            
            Returns:
            --------
            list HDUs
            '''
                pass

            def _remove_scattered_light_problem_areas(self,header , imageData, **kwargs):
            '''
                TESS has certain regions (of the chip(s)) in certain exposures that are known to have ...
                ... high levels of polluting scattered light
                We may want to completely mask these out
                
                Input:
                --------
                list HDUs
                
                Returns:
                --------
                list HDUs
                '''
            pass

    def _remove_strap_regions(self,header , imageData, **kwargs):
        '''
            TESS has certain regions (of the chip(s)) in which the backing-straps provide confusingly high signals
            We may want to completely mask these out (or perhaps flag them in some alternative manner)
            
            Input:
            --------
            list HDUs
            
            Returns:
            --------
            list HDUs
            '''
                pass








