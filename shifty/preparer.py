'''
   Classes / methods to prepare fits files
   Provides methods to 
   (i) "clean" data and 
   (ii) combine multiple image-sets into a single, large, "stacked" fits-file
'''

# -------------------------------------------------------------------------------------
# Third party imports
# -------------------------------------------------------------------------------------
import os, sys
from collections import OrderedDict
import numpy as np
import copy

import astropy
from astropy.io import fits
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

class ImagePreparer(Downloader):
    '''
        (1)Loads fits-files
        (2)Cleans/Prepares data from fits files
        (3)Instantiates "stack" fits-file
        
        Parent class for ...
        - TessImagePreparer, HubbleImagePreparer, PanstarrsImagePreparer, ...
        
        inputs:
        -------
        None
        
        methods:
        --------
        _load_image()
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
    def generate_cleaned_stack_file(self,):
        '''
            This function will:
            - Get data from file(s)
            - Do "cleaning"
            - Save all component data into a single large fits file
            
            *** STUB FUNCTION THAT WILL BE OVERWRITTEN BY CHILD CLASS ***
            
            Input:
            ------
            
            Returns:
            --------
            '''
        pass
    
    
    
    # -------------------------------------------------------------------------------------
    # The methods below are for the loading of *general* fits-format data files
    # -------------------------------------------------------------------------------------
    
    def _load_image(self , fits_filepath):
        '''
            Load a single image
            - Currently a wrapper around astropy.fits.open
            - With the potential for added functionality (to be added later)
            
            Input:
            ------
            fits_filepath
             - valid filepath to single, valid fits-file
            
            Returns:
            --------
            lastropy.io.fits.hdu.hdulist.HDUList object
            - [[ currently defaults to "None" if filepath invalid]]
            
            '''
        # Here I am returning the HDUlist object
        # ** IT NEEDS TO BE CLOSED LATER ON !! **
        return fits.open(fits_filepath) if os.path.isfile(fits_filepath) and '.fits' in fits_filepath else None
    
    






class TESSImagePreparer(ImagePreparer, TESSDownloader):
    '''
        For the specific preparation of *TESS* data 
        
        (1)Loads fits-files
        (2)Cleans/Prepares data from fits files
        (3)Instantiates "stack" fits-file
        
        inputs:
        -------
        None
        
        methods:
        --------
        _remove_stars()  (i.e. set pixels that contain stars to NaN)
        _remove_bad_cadences()
        _remove_scattered_light_problem_areas()
        _remove_strap_regions()
        
        main public method:
        -------------------
        get_image_data_set() => Returns an ImageDataSet
        
    '''
    
    def __init__(self, ) :
        
        # - Allow ourselves to inherit methods
        ImagePreparer.__init__(self, )
        TESSDownloader.__init__(self, )
        
        # - Define some important TESS-related quantities
        self.obs_code = 'C57'
    
    # -------------------------------------------------------------------------------------
    # Public Methods
    # -------------------------------------------------------------------------------------
    
    def generate_cleaned_data(self, **kwargs ):
        '''
            This function will:
             - read raw component fits files
             - allow all aspects of "cleaning"
             - save all component data into a ...
            
            returns:
            --------
            ...
        '''
        
        # Parse the file-spec and decide what fits-files will be loaded
        fits_filepaths = self._parse_filespec( **kwargs )
        x0,x1,y0,y1    = self._parse_patchspec( **kwargs )
        
        # Initialize a master HUDlist for a new "stacked" fits-file
        # - Includes nicely filled-out primary-header
        stack_fits_filepath = os.path.join(self._fetch_data_directory(), 'shift.fits')
        new_hdul = self._initialize_stack_HDUlist(stack_fits_filepath, **kwargs)

        # Loop over the files sequentially
        for fits_filepath in fits_filepaths:
            
            # open & read the individual TESS fits file
            with fits.open(fits_filepath) as hdul:
                header, imagedata, unc = hdul[1].header , hdul[1].data, hdul[2].data

                # clean the data
                clean_imagedata = self._clean_data(header, imagedata, **kwargs )
                
                # create an ImageHDU object
                # make ImageHDU have header & data from the TESS fits file we are reading
                # add to ImageHDU stack-file
                # [[ Note that I am not doing anything with the uncertainty data at present ]]
                new_hdul.append( fits.ImageHDU( header=hdul[1].header , data=clean_imagedata ) )

        # save it all to stacked fits-file
        new_hdul.writeto(stack_fits_filepath, overwrite=True)
        
        # ensure that the new stack-file is closed ...
        new_hdul.close()
        del new_hdul
        
        return stack_fits_filepath



    def get_prfs(self, ):
        '''
            ...
        '''
        # Defined in downloader ...
        # - N.B.: if the prfs already exist, *get_prf()* does no download work ...
        # ... but instead returns list of available filepaths
        prf_filepaths = self.download_prf()
    
        # Use standard method to open filepaths and return list of HDUs
        #return self._load_images(fits_filepaths = prf_filepaths)
    
    
    
    # -------------------------------------------------------------------------------------
    # The method(s) below are for the deciding which TESS fits-format data files to load
    # -------------------------------------------------------------------------------------
    
    def _parse_filespec(self,  **kwargs ):
        '''
            parsing function to allow passing of a variety of arguments to ...
            ... the get_image_data_set() function
            
            Thus far it knows how to interpret
            (i) a list of fits-filepaths
            (ii) a request for 'development' fits-files
            (iii) a request for a specific sector/camera/chip
            
            '''
        
        try:
            
            # if explict filepaths defined, use this, and pass on
            if 'fits_filepaths' in kwargs:
                fits_filepaths = [ ffp for ffp in np.atleast_1d(kwargs['fits_filepaths']) if os.path.isfile(ffp) and '.fits' in ffp ]
            
            # If 'development' is specified, then get limited test/development data-set
            elif  'development' in  kwargs and kwargs['development']:
                fits_filepaths = self._ensure_test_data_available_locally()
            
            # If the sector/camera/chip specified, then get the required filepaths
            elif    np.all( [_ in kwargs for _ in ['sectorNumber', 'cameraNumber', 'chipNumber']] ) \
                and isinstance(kwargs['sectorNumber'], int) \
                    and isinstance(kwargs['cameraNumber'], int) \
                        and isinstance(kwargs['chipNumber'], int):
                directory_path = os.path.join( self.tess_dir,
                                              str(kwargs['sectorNumber']),
                                              str(kwargs['cameraNumber']),
                                              str(kwargs['chipNumber']))
                fits_filepaths = glob.glob( os.path.join(directory_path , '*.fits') )
            else:
                pass
    
        except Exception as error:
            print('Could not interpret the supplied argument to get_image_data_set() : see _parse_filespec()')
            print('*** NO FILES WILL BE OPENED *** ')
            print(error)
            fits_filepaths = []

        return fits_filepaths
            
            
    def _parse_patchspec(self, **kwargs):
        '''
            parsing function to allow passing of a variety of means to ...
            ... specify the sub-region (patch) of a chip to work with
            
            Thus far it knows how to interpret
            (i) pythonic, 0-based, array specification
            (ii) pixel, 1-based, specification
        '''
        x0,x1,y0,y1 = 0,-1,0,-1
        try:
            if 'patch' in kwargs and kwargs['patch']:
        
                # If python-like specified, then array-elements numbered from 0
                if  'python' in kwargs and kwargs['python'] == True and 'xlim' in kwargs and 'ylim' in kwargs:
                    x0,x1,y0,y1 = kwargs['xlim'][0], kwargs['xlim'][1], kwargs['ylim'][0], kwargs['ylim'][1]
    
                # If fits-pixel-like specified, then need to offset
                elif 'pixel' in kwargs and kwargs['pixel'] == True and 'xlim' in kwargs and 'ylim' in kwargs:
                    x0,x1,y0,y1 = kwargs['xlim'][0]-1, kwargs['xlim'][1]-1, kwargs['ylim'][0]-1, kwargs['ylim'][1]-1
    
                else:
                    pass
    
        except Exception as error:
            print('Could not parse patch-specification')
            print(error)

        return x0,x1,y0,y1

    # -------------------------------------------------------------------------------------
    # The method(s) below are convenience functions while developing ...
    # -------------------------------------------------------------------------------------
  
    def _load_test_images(self,):
        ''' 
            Convenience function to load a small, pre-defined sample of test data
            Does *NOT* create a stack file, just opens a bunch of individual fits-files
        '''
        return [ self._load_image(fp) for fp in self._ensure_test_data_available_locally() ]

    def _generate_test_stack_file(self,):
        '''
            Convenience function to load a small, pre-defined sample of test data into a single stack-file
            Does *NOT* do any cleaning
        '''
        return self.generate_cleaned_stack_file( development = True )

    """
    def _load_sector_camera_chip(self, sectorNumber, cameraNumber, chipNumber):
    ''' convenience function to load data for single sector/camera/chip'''
    return self._load_images( { 'sectorNumber' : sectorNumber,
    'cameraNumber' : cameraNumber,
    'chipNumber'   : chipNumber} )
    """
    
    # -------------------------------------------------------------------------------------
    # The method(s) below are for creating/handling an overall "stack" fits-file
    # -------------------------------------------------------------------------------------
    
    def _initialize_stack_HDUlist(self, stack_filepath, **kwargs ):
        ''' 
            This will hold the 'stacked' fits-file
            
            [[MIGHT WANT TO RECORD DETAILS OF **kwargs INTO PRIMARY HEADER: E.g. RECORD CLEANING METHOD, ETC]]
            
            [[Might want to move this to the parent class]]
        '''
        # Create the HDUlist
        stack_hdulist = fits.HDUList([fits.PrimaryHDU()])
        
        # Populate the header [[not yet done ]]
        
        # Save to file and return HDUlist
        stack_hdulist.writeto(stack_filepath, overwrite=True)
        return stack_hdulist


    # -------------------------------------------------------------------------------------
    # The methods below are for the "PARSING" of TESS HDUS from FITS
    # -------------------------------------------------------------------------------------

    
    def _get_midtime(self, header ):
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
        
        try:
            # Check time system is Barycentric Dynamical Time (TDB)
            assert 'TDB' in header['TIMESYS'], \
                'wrong TIMESYS ... '
            
            # Calculate the exposure mid-time
            T = Time(       header['BJDREFI']    \
                     +      header['BJDREFF']    \
                     + 0.5*(header['TSTART']      + header['TSTOP']) ,
                     format='jd', scale='tdb')
    
        except Exception as error :
            print('There was an error calculating the exposure mid-time')
            print(error)
            T = None
        
        return T
    
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
    
        return imageData



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



class HSTImagePreparer(ImagePreparer):
    '''
        You know we'll want to do it at some point !!!
        '''
    pass








