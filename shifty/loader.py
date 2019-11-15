'''
   Classes / methods to load fits file from disk into an "ImageDataSet" object
'''

# -------------------------------------------------------------------------------------
# Third party imports
# -------------------------------------------------------------------------------------
import os, sys
import astropy
from astropy.io import fits
from collections import OrderedDict
import numpy as np
from astropy.time import Time

# -------------------------------------------------------------------------------------
# Any local imports
# -------------------------------------------------------------------------------------
from downloader import *
import data
from data import ImageDataSet
from refcat import RefCat

# -------------------------------------------------------------------------------------
# Various class definitions for *data import * in shifty
# -------------------------------------------------------------------------------------

class ImageLoader(Downloader):
    '''
        (1)Loads fits-files
        (2)Cleans/Prepares data from fits files
        (3)Instantiates ImageDataSet object
        
        Parent class for ...
        - TessImageLoader, HubbleImageLoader, PanstarrsImageLoader, ...
        
        inputs:
        -------
        ???sector number, camera, chip, detector_range, bad_data_thresholds ???
        
        methods:
        --------
        _load_images()
        _remove_stars()  (i.e. set pixels that contain stars to NaN)
        _remove_bad_cadences()
        _remove_scattered_light_problem_areas()
        _remove_strap_regions()
        
        main public method:
        -------------------
        get_image_data_set() => Returns an ImageDataSet
        
        
        '''
    
    def __init__(self, ) :
        
        # - Allow ourselves to use Downloader methods
        Downloader.__init__(self, )
        
        # - Local directory for saving data
        self.local_dir = self._fetch_data_directory()
    
    # -------------------------------------------------------------------------------------
    # Public Methods
    # -------------------------------------------------------------------------------------
    def get_image_data_set(self,):
        '''
            Overall image loader as envisaged by Geert & Matt over coffee
            - Gets data from file(s)
            - Does "cleaning"
            - Creates ImageDataSet object
            
            *** STUB FUNCTION THAT WILL BE OVERWRITTEN BY CHILD CLASS ***
            
            Input:
            ------
            (1) ???? list of valid filepaths to valid fits-files ????
            (2) option args to pass through to various sub-functions
            
            Returns:
            --------
            ImageDataSet
            '''
        pass
    
    
    
    # -------------------------------------------------------------------------------------
    # The methods below are for the loading of *general* fits-format data files
    # -------------------------------------------------------------------------------------
    
    def _load_images(self, *args, **kwargs ):
        '''
            Load multiple images
            (1) Interprets file_spec_container to decide what files need to be loaded
            - ADDITIONAL PRE-FILTERING CAN/WILL BE DONE BY CHILD METHODS
            (2) Uses "_load_image()" to open the files
            
            Input:
            ------
            valid filepath(s) to valid fits-files
            
            Returns:
            --------
            list of astropy.io.fits.hdu.hdulist.HDUList objects
            '''
        
        # if kwargs contains "fits_filepaths" that are valid, then we are good to go ...
        if 'fits_filepaths' in kwargs:
            try:
                fits_filepaths = [ ffp for ffp in np.atleast_1d(kwargs['fits_filepaths']) if os.path.isfile(ffp) and '.fits' in ffp ]
            except Exception as error:
                fits_filepaths = []
                print('problem parsing fits_filepaths : %r' % kwargs['fits_filepaths'] )
                print(error)
    
        # no other handling-methods currently in place (but CHILD may have PRE-FILTERED)
        else:
            print(' *** At present, no method is in place to interpret this input *** ')
            print(' ******              No files will be loaded                 ******')
            fits_filepaths = []
        
        # open filepaths
        return [ self._load_image(fp) for fp in fits_filepaths ]
    
    def _load_image(self , fits_filepath):
        '''
            Load a single image
            - Currently a wrapper around astropy.fits.open
            - With the potential for added functionality (to be added later)
            
            Input:
            ------
            valid filepath to valid fits-file
            
            Returns:
            --------
            lastropy.io.fits.hdu.hdulist.HDUList object
            - [[ currently defaults to "None" if filepath invalid]]
            
            '''
        return fits.open(fits_filepath) if os.path.isfile(fits_filepath) and '.fits' in fits_filepath else None








class TESSImageLoader(ImageLoader, TESSDownloader):
    '''
        
        => input: sector number, camera, chip, detector_range, bad_data_thresholds
        
        _load_images()
        _remove_stars()  (i.e. set pixels that contain stars to NaN)
        _remove_bad_cadences()
        _remove_scattered_light_problem_areas()
        _remove_strap_regions()
        get_image_data_set() => ImageDataSet
        
        '''
    
    def __init__(self, ) :
        
        # - Allow ourselves to inherit methods
        ImageLoader.__init__(self, )
        TESSDownloader.__init__(self, )
        
        # - Define some important TESS-related quantities
        self.obs_code               = 'C57'
    
    # -------------------------------------------------------------------------------------
    # Public Methods
    # -------------------------------------------------------------------------------------
    def get_image_data_set(self, **kwargs ):
        '''
            Overall image loader
            - Gets data from file(s)
            - Does "cleaning"
            - Creates ImageDataSet object
            
            Inputs:
            ------
            (1) params to pass through to _load_images
            (2) params to pass through to _clean_data
            
            Returns:
            --------
            ImageDataSet
            '''
        # Load the data from files
        HDUs = self._load_images( **kwargs )
            
        # parse the HDUs to get requisite data to create ImageDataSet object
        # [[may be more useful to move this up before the cleaning stage]]
        HDU_wcs_headers, HDU_data, HDU_unc, HDU_midtimes = self._parse_HDUs_for_ImageDataSet( HDUs )

        # clean the data
        # - [[could just ovewrite original HDUs with clean version: will leave both for now]]
        clean_HDU_data = self._clean_data(HDU_wcs_headers, HDU_data, **kwargs )

        # create & return ImageDataSet
        return ImageDataSet(HDU_wcs_headers, clean_HDU_data, HDU_midtimes, self.obs_code, HDU_unc = HDU_unc)
    
    
    def get_prfs(self, ):
        '''
            ...
        '''
        # Defined in downloader ...
        # - slight hack: down_prf returns list of downloaded filepaths
        #  (but does no download work if the prfs already exist)
        prf_filepaths = self.download_prf()
    
        # Use standard method to open filepaths and return list of HDUs
        return self._load_images(fits_filepaths = prf_filepaths)
    
    # -------------------------------------------------------------------------------------
    # The methods below are for the "PARSING" of TESS HDUS from FITS
    # -------------------------------------------------------------------------------------
    def _parse_HDUs_for_ImageDataSet(self, HDUs):
        '''
            Parse a set of TESS HDUs 
             - splitting out into lists / arrays of ...
             (1) header data containing the WCS
             (2) the image data 
             (3) the exposure mid-times
        '''
        HDU_wcs_headers = [_[1].header for _ in HDUs]
        HDU_data        = np.array( [_[1].data for _ in HDUs] )
        HDU_unc         = np.array( [_[2].data for _ in HDUs] )
        HDU_midtimes    = self._get_midtimes( HDU_wcs_headers )
        return HDU_wcs_headers, HDU_data, HDU_unc, HDU_midtimes,
    
    def _get_midtimes(self, header_iterable ):
        '''
            Get a BJD/TDB out of the headers for each exposure
            
            Input:
            ------
            list/iterable of fits.header objects
            
            Returns:
            --------
            astropy.time.core.Time object
        '''
        # (1) check that it is always TDB in the 'cal' header [TIMESYS]
        # (2) LIVETIME gets effective exposure time [[not used currently]]
        # (3) BARYCORR might provide a useful sanity check
        # N.B.
        # ['BJDREFI'] / integer part of BTJD reference date
        # ['BJDREFF'] / fraction of the day in BTJD reference date
        # ['TSTART']  / observation start time in BTJD
        # ['TSTOP']   / observation stop time in BTJD
        t_arr = []
        
        try:
            for h in header_iterable:
                # Check time system is Barycentric Dynamical Time (TDB)
                assert 'TDB' in h['TIMESYS'], 'wrong TIMESYS ... '
                # Calculate the exposure mid-time
                t_arr.append(       h['BJDREFI']    \
                             +      h['BJDREFF']    \
                             + 0.5*(h['TSTART']      \
                                    +      h['TSTOP'])    )
            
            T = Time(t_arr, format='jd', scale='tdb')
        except Exception as error :
            print('There was an error calculating the exposure mid-times')
            print(error)
            T = None
        
        return T
    
    # -------------------------------------------------------------------------------------
    # The methods below are for the "CLEANING" of TESS data
    # -------------------------------------------------------------------------------------
    def _clean_data(self, HDU_wcs_headers, HDU_data, **kwargs):
        '''
            Wrapper around the various "cleaning" methods below
            - Simply provided as a means to enable simple toggling on/off of functions
            - Uses cleaning_parameters to control what is/isn't evaluated
            '''
        # dict to hold key:function mappings
        # - Have used ORDERED because at some point it might be important what order the funcs are used
        cleaning_function_dict = OrderedDict([
                                              ('mask'      , self._mask_stars ),
                                              ('subtract'  , self._subtract_stars ),
                                              ('bad_cad'   , self._remove_bad_cadences ),
                                              ('scat'      , self._remove_scattered_light_problem_areas ),
                                              ('strap'     , self._remove_strap_regions ),
                                              ])
            
        # loop over possible funcs (in order of dict)
        # [[ NOTICE THE IMPLICIT DESIGN CHOICE THAT ALL CLEANING FUNCTIONS MUST RETURN HDU_data]]
        for key, func_to_run in cleaning_function_dict.items():
            # run a function if it is included as True in cleaning_parameters (e.g. {'mask':True})
            if key in kwargs and kwargs[key]:
                HDU_data = func_to_run(HDU_wcs_headers , HDU_data, **kwargs)
    
        return HDU_data



    def _mask_stars(self, HDU_wcs_headers, HDU_data, **kwargs):
        '''
            We want to remove stars in some way
            Barentsen & Payne discussed a simple mask: i.e. set pixels that contain stars to NaN
            This would be done based on GAIA positions
            
            This is *NOT* subtraction (see _subtract_stars below )
            
            Presumably only one of _mask_stars / _subtract_stars is required, but I am 100% certain that Holman will at least want to experiment with subtraction
            
        '''
        
        # provide a means to only do the refcat search once
        # - using this assumes that all of the images are closely aligned (v. similar ra,dec ranges)
        ONESHOT = True if 'oneshot' in kwargs and kwargs['oneshot'] else False
        
        if ONESHOT:
            # find the location of all of the stars on the FIRST image ONLY
            ra, dec , pix = RefCat().find_all_stars_on_image(HDU_wcs_headers[0], HDU_data[0])
        
        for header, image_data in zip(HDU_wcs_headers, HDU_data):
            
            if not ONESHOT:
                # find the location of all of the stars on each individual image
                ra, dec , pix = RefCat().find_all_stars_on_image(header, image_data)
     
            # presumably need to do something about deciding how big a mask to use
            # - based on the source magnitude?
            # Perhaps something from photutils
            # https://photutils.readthedocs.io/en/stable/psf.html
            # http://docs.astropy.org/en/stable/api/astropy.convolution.discretize_model.html
            print(' ** WARNING: just outputing a single central mask pixel at present ** ')
            
            # mask all of the stars
            image_data[pix] = 0
            
        return HDU_data

    def _subtract_stars(self,HDUs, **kwargs):
        '''
            We want to remove stars in some way
            Holman & Payne have generally assumed some form of subtraction
            
            This is *NOT* masking (see _mask_stars )
            
            Presumably only one of _mask_stars / _subtract_stars is required, but I am 100% certain that Holman will at least want to experiment with subtraction
            
            Input:
            --------
            list HDUs
            
            Returns:
            --------
            list HDUs
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
        
        
        return HDUs

    def _remove_bad_cadences(self,HDUs, **kwargs):
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
        return HDUs

    def _remove_scattered_light_problem_areas(self,HDUs, **kwargs):
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
        return HDUs

    def _remove_strap_regions(self,HDUs, **kwargs):
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
        return HDUs



    # -------------------------------------------------------------------------------------
    # The method(s) below are for the loading of TESS fits-format data files
    # -------------------------------------------------------------------------------------
    
    def _load_images(self, *args, **kwargs):
        '''
             This child method does PRE-FILTERING before handing off to ...
             ... PARENT *_load_images()* method in ImageLoader
             - PARENT knowns how to handle filepaths to fits_files
             
             Input:
             ------
             valid file_spec_container
             
             Returns:
             --------
             list of astropy.io.fits.hdu.hdulist.HDUList objects
        '''
        to_parent = {}

        # if explict filepaths defined, use this, and pass through to parent
        if 'fits_filepaths' in kwargs:
            try:
                to_parent['fits_filepaths'] = [ ffp for ffp in np.atleast_1d(kwargs['fits_filepaths']) if os.path.isfile(ffp) and '.fits' in ffp ]
            except Exception as error:
                to_parent['fits_filepaths'] = []

        # If 'development' is specified, then get limited test/development data-set
        elif  'development' in  kwargs:
            # Use the test-data method to ensure files available locally
            to_parent['fits_filepaths'] = self._ensure_test_data_available_locally()

        # If the sector/camera/chip specified, then get the required filepaths
        # [[ *** AS WRITTEN THIS WILL ONLY RETURN THE DATA FOR A SINGLE CHIP *** ]]
        elif    np.all( [_ in kwargs for _ in ['sectorNumber', 'cameraNumber', 'chipNumber']] ) \
                and isinstance(kwargs['sectorNumber'], int) \
                and isinstance(kwargs['cameraNumber'], int) \
                and isinstance(kwargs['chipNumber'], int):
            directory_path = os.path.join( self.tess_dir,str(sectorNumber),str(cameraNumber),str(chipNumber))
            try:
                to_parent['fits_filepaths'] = glob.glob( os.path.join(directory_path , '*.fits') )
            except Exception as error:
                to_parent['fits_filepaths'] = []

        # Might be something that the parent knows how to parse ...
        else:
            to_parent = kwargs

        return super(TESSImageLoader, self)._load_images( **to_parent  )


    def _load_test_images(self,):
        ''' convenience function to load a small, pre-defined sample of test data'''
        return self._load_images( development = True )
    
    def _load_sector_camera_chip(self, sectorNumber, cameraNumber, chipNumber):
        ''' convenience function to load data for single sector/camera/chip'''
        return self._load_images( { 'sectorNumber' : sectorNumber,
                                    'cameraNumber' : cameraNumber,
                                    'chipNumber'   : chipNumber} )


class HSTImageLoader(ImageLoader):
    '''
        You know we'll want to do it at some point !!!
        '''
    pass








