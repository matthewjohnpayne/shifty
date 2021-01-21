# -*- coding: utf-8 -*-
# shifty/shifty/imagehandler.py

'''
   Classes / methods for reading cleaned image fits files.
   Provides methods to
   -
'''

# -----------------------------------------------------------------------------
# Third party imports
# -----------------------------------------------------------------------------
import os
import sys
from collections import OrderedDict
import numpy as np

from astropy.io import fits
from astropy.time import Time

# -----------------------------------------------------------------------------
# Any local imports
# -----------------------------------------------------------------------------
sys.path.append(os.path.dirname(os.path.dirname(
                os.path.realpath(__file__))))
from shifty.downloader import Downloader
from shifty.imagedata import ImageDataSet
from shifty.refcat import RefCat

# -----------------------------------------------------------------------------
# Various class definitions for *data import * in shifty
# -----------------------------------------------------------------------------

class ImageReader(Downloader):
    '''
        (1)Loads fits-files
        (2)Reads WCS
        (3)Calculate pixel coordinates
        (4)Interact with 'known'?

        inputs:
        -------
        Needs to somehow get knowledge of which Header Keywords to look for

        methods:
        --------
        loadImageAndHeader()
        _find_key_value()

        main public method:
        -------------------
        loadImageAndHeader()


        '''

    def __init__(self, **kwargs):
        # - Allow ourselves to use Downloader methods
        Downloader.__init__(self,)
        # - Local directory for saving data
        self.local_dir = self._fetch_data_directory()

    # -------------------------------------------------------------------------
    # Public Methods
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # The methods below are for the loading of *general* fits-format data files
    # -------------------------------------------------------------------------

    def loadImageAndHeader(self, inputFile, extno=None, verbose=False, **kwargs):
        '''
        Reads in a fits file (or a given extension of one).
        Returns the image data, the header, and a few useful keyword values.

        input:
        ------
        fits_filepath  - str      - valid filepath to single, valid fits-file
        extno          - int OR None  - Extension number that has the image data (if multi-extension)
        verbose        - bool         - Print extra stuff if True
        EXPTIME        - str OR float - Exposure time in seconds (keyword or value)
        MAGZERO        - str OR float - Zeropoint magnitude (keyword or value)
        MJD_START      - str OR float - MJD at start of exposure (keyword or value)
        GAIN           - str OR float - Gain value (keyword or value)
        FILTER         - str          - Filter name (keyword)
        NAXIS1         - str OR int   - Number of pixels along axis 1
        NAXIS2         - str OR int   - Number of pixels along axis 2
        INSTRUMENT     - str          - Instrument name (keyword)
        
        A float value can be defined for most keywords, rather than a keyword name;
        this will use that walue rather than searching for the keyword in the headers.  
        INSTRUMENT and FILTER obviously can't be floats,
        so use a leading '-' to specify a value rather than a keyword name to use. 

        output:
        -------
        self.header_keywords - dictionary - A bunch of header keyword names. 
                                          - Not sure these are needed after use here.
        self.EXPTIME         - float      - Exposure time in seconds
        self.MAGZERO         - float      - Zeropoint magnitude
        self.MJD_START       - float      - MJD at start of exposure
        self.MJD_MID         - float      - MJD at centre of exposure
        self.GAIN            - float      - Gain value
        self.FILTER          - str        - Filter name
        self.NAXIS1          - int        - Number of pixels along axis 1
        self.NAXIS2          - int        - Number of pixels along axis 2
        '''
        # Initialization of default standard Fits header keywords
        # - These may be overwritten by kwargs
        self.header_keywords = {
                                'EXPTIME': 'EXPTIME',      # Exposure time in seconds
                                'MAGZERO': 'MAGZERO',      # Zeropoint magnitude
                                'MJD_START': 'MJD-OBS',    # MJD at start of exposure
                                'GAIN': 'GAIN',            # Gain value
                                'FILTER': 'FILTER',        # Filter name
                                'NAXIS1': 'NAXIS1',        # Number of pixels along axis 1
                                'NAXIS2': 'NAXIS2',        # Number of pixels along axis 2
                                'INSTRUMENT': 'INSTRUME',  # Instrument name
                                }
        # This is not a comprehensive list, just the ones that I expect we'll need,
        # particularly those that I know have different names from different telescopes. 

        # Do a loop over the kwargs and see if any header keywords need updating
        # (becasue they were supplied)
        for key, value in kwargs.items():
            if key in self.header_keywords:
                self.header_keywords[key] = value

        # Read the file:
        with pyf.open(inputFile) as han:
          if extno is None:    # Do this if it's a single extension file (single chips)
            print('Warning: Treating this as a single extension file.')
            self.data = han.data
            self.header = han.header
            self.header0 = self.header   # For single extension, header and header0 are the same.
          else:                # Do this for multi-extension files (mosaic cameras)
            self.data = han[extno].data
            self.header = han[extno].header  # Header for the particular extension
            self.header0 = han[0].header  # Extension 0 usually holds an general header for mosaic

        # Read the WCS
        self.WCS = wcs.WCS(self.header)

        # Use the defined keywords to save the values into self.
        # Search both headers if neccessary (using keyValue)
        for key, use in self.header_keywords:
            if key not in ['INSTRUMENT', 'FILTER', 'NAXIS1', 'NAXIS2']:
                # Most keywords can just have a float value defined instead of keyword name
                setattr(self, key, float(self._find_key_value(use)) if type(use) != float else use)
            elif key in ['NAXIS1', 'NAXIS2']:
                # Some keywords can have an integer value defined instead of keyword name
                setattr(self, key, int(self._find_key_value(use)) if type(use) != int else use)
            elif key == 'INSTRUMENT':
                # INSTRUMENT and FILTER obviously can't be floats,
                # so use a leading '-' to specify a value rather than a keyword name to use. 
                self.INSTRUMENT = self._find_key_value(use)) if use[0] != '-' else use[1:]
            elif key == 'FILTER':
                # Filter only uses the first character of the supplied, not all the junk
                self.FILTER = self._find_key_value(use))[0] if use[0] != '-' else use[1]
        
        # Also define the middle of the exposure:
        self.MJD_MID = self.MJD_START + self.EXPTIME / 172800.0

        print('{}\n'.format((self.EXPTIME, self.MAGZERO, self.MJD_START, self.MJD_MID,
                             self.GAIN, self.FILTER, self.NAXIS1, self.NAXIS2, 
                             self.INSTRUMENT)) if verbose else '', end='')

    def _find_key_value(self, key):
        """
        First checks extension header for a keyword; if fails, checks main header.
        This is neccessary because, annoyingly, some telescopes put things like the EXPTIME 
        in the main header, while others put it in the header for each extension and NOT the main one.
        
        input:
        ------
        key - str - keyword to look for

        output:
        -------
        value of keyword found in headers
        """
        try:
          value =self.header[key]
        except KeyError:
          value = self.header0[key]
        return value


    def generate_cleaned_data(self, **kwargs):
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
        fits_filepaths = self._parse_filespec( **kwargs)
        x0,x1,y0,y1    = self._parse_patchspec( **kwargs)

        # Initialize a master HUDlist for a new "stacked" fits-file
        # - Includes nicely filled-out primary-header
        stack_fits_filepath = os.path.join(self._fetch_data_directory(), 'shift.fits')
        new_hdul = self._initialize_stack_HDUlist(stack_fits_filepath, **kwargs)

        # Loop over the files sequentially
        for fits_filepath in fits_filepaths:

            # open & read the individual TESS fits file
            with fits.open(fits_filepath) as hdul:
                header, imagedata, unc = hdul[1].header, hdul[1].data, hdul[2].data

                # clean the data
                clean_imagedata = self._clean_data(header, imagedata, **kwargs)

                # create an ImageHDU object
                # make ImageHDU have header & data from the TESS fits file we are reading
                # add to ImageHDU stack-file
                # [[ Note that I am not doing anything with the uncertainty data at present ]]
                new_hdul.append( fits.ImageHDU( header=hdul[1].header, data=clean_imagedata))

        # save it all to stacked fits-file
        new_hdul.writeto(stack_fits_filepath, overwrite=True)

        # ensure that the new stack-file is closed ...
        new_hdul.close()
        del new_hdul

        return stack_fits_filepath


    # -------------------------------------------------------------------------
    # The method(s) below are convenience functions while developing ...
    # -------------------------------------------------------------------------



    # -------------------------------------------------------------------------
    # The method(s) below are for creating/handling an overall "stack" fits-fil
    # -------------------------------------------------------------------------

    def _initialize_stack_HDUlist(self, stack_filepath, **kwargs):
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


    # -------------------------------------------------------------------------
    # The methods below are for the "PARSING" of TESS HDUS from FITS
    # -------------------------------------------------------------------------


    def _get_midtime(self, header):
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
                     + 0.5*(header['TSTART']      + header['TSTOP']),
                     format='jd', scale='tdb')

        except Exception as error:
            print('There was an error calculating the exposure mid-time')
            print(error)
            T = None

        return T

    # -------------------------------------------------------------------------
    # The methods below are for the "CLEANING" of TESS data
    # -------------------------------------------------------------------------
    def _clean_data(self, header, imageData, **kwargs):
        '''
            Wrapper around the various "cleaning" methods below
            - Simply provided as a means to enable simple toggling on/off of functions
            - Uses kwargs to control what is/isn't evaluated
            '''
        # dict to hold key:function mappings
        # - Have used ORDERED because at some point it might be important what order the funcs are used
        cleaning_function_dict = OrderedDict([
                                              ('mask'     , self._mask_stars),
                                              ('subtract' , self._subtract_stars),
                                              ('clip'     , self._clip_peaks),
                                              ('bad_cad'  , self._remove_bad_cadences),
                                              ('scat'     , self._remove_scattered_light_problem_areas),
                                              ('strap'    , self._remove_strap_regions),
                                              ])

        # loop over possible funcs (in order of dict)
        # [[ N.B. HDU_data should be modified in-place by functions ]]
        for key, func_to_run in cleaning_function_dict.items():
            # run a function if it is included as True in cleaning_parameters (e.g. {'mask':True})
            if key in kwargs and kwargs[key]:
                imageData = func_to_run(header, imageData, **kwargs)

        return imageData



    def _mask_stars(self, header, imageData, **kwargs):
        '''
            We want to remove stars in some way
            Barentsen & Payne discussed a simple mask: i.e. set pixels that contain stars to NaN
            This would be done based on GAIA positions

            This is *NOT* subtraction (see _subtract_stars below)

            Presumably only one of _mask_stars / _subtract_stars / _clip_peaks is required

        '''

        # Provide a means to only do the refcat search once
        # If NO useful refcat dictionary supplied, do search, otherwise use supplied dictionary
        # - using this implicitly assumes that all of the images are closely aligned (v. similar ra,dec ranges)
        # -
        if 'refcat_dict' not in kwargs or kwargs['refcat_dict'] == {}:
            print('in calc loop')
            kwargs['refcat_dict'] = {}
            ra,dec,pix,int_pix = RefCat().find_all_stars_on_image(header, imageData)
            kwargs['refcat_dict']['ra'], kwargs['refcat_dict']['dec'], kwargs['refcat_dict']['pix'], kwargs['refcat_dict']['int_pix'] = ra,dec,pix, int_pix


        # Need to do something about deciding how big a mask to use, based on the source magnitude
        # Perhaps something from photutils
        # https://photutils.readthedocs.io/en/stable/psf.html
        # http://docs.astropy.org/en/stable/api/astropy.convolution.discretize_model.html
        # Perhaps using the downloaded prf
        print(' ** WARNING: just outputting a single central mask pixel at present ** ')


        # mask all of the stars
        # - N.B. this is likely to fail if nPixels > 0 in RefCat().find_all_stars_on_image()
        # - N.B. this alters imageData-in-place ...
        rows, cols = kwargs['refcat_dict']['int_pix'][1], kwargs['refcat_dict']['int_pix'][0]
        imageData[rows, cols] = 0


    def _subtract_stars(self, header, imageData, **kwargs):
        '''
            We want to remove stars in some way
            Holman & Payne have generally assumed some form of subtraction

            This is *NOT* masking (see _mask_stars)

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

    def _clip_peaks(self, header, imageData, **kwargs):
        '''
            We want to remove stars in some way
            Could use some form of peak-clipping to remove any/all points above some ~background level

            Presumably only one of _mask_stars / _subtract_stars / _clip_peaks is required

            Input:
            --------

            Returns:
            --------
            '''


    def _remove_bad_cadences(self,header, imageData, **kwargs):
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

    def _remove_scattered_light_problem_areas(self,header, imageData, **kwargs):
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

    def _remove_strap_regions(self,header, imageData, **kwargs):
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
