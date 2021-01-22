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
import numpy as np

from astropy.io import fits
from astropy.time import Time
from astropy import wcs

# -----------------------------------------------------------------------------
# Any local imports
# -----------------------------------------------------------------------------
sys.path.append(os.path.dirname(os.path.dirname(
                os.path.realpath(__file__))))
from shifty.downloader import Downloader

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
    # Turning off stupid some syntax-checker warnings:
    # pylint: disable=too-many-instance-attributes
    # Why on earth should an object only have 7 attributes?
    # That seems dumb. Turning off this warning.
    # pylint: disable=attribute-defined-outside-init
    # Why on earth should attributes not be defined outside init?
    # That seems dumb. Turning off this warning.

    def __init__(self, filename=None, extno=None, verbose=False, **kwargs):
        # - Allow ourselves to use Downloader methods
        Downloader.__init__(self,)
        # - Local directory for saving data
        self.local_dir = self._fetch_data_directory()
        self.header_keywords = {}
        self.keyword_values = {}
        self.filename = filename
        self.extno = extno
        if filename:
            self.loadImageAndHeader(filename, extno, verbose=verbose, **kwargs)

    # -------------------------------------------------------------------------
    # Public Methods
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # The methods below are for the loading of *general* fits-format data files
    # -------------------------------------------------------------------------

    def loadImageAndHeader(self, filename=None, extno=None, verbose=False,
                           **kwargs):
        '''
        Reads in a fits file (or a given extension of one).
        Returns the image data, the header, and a few useful keyword values.

        input:
        ------
        filename       - str          - valid filepath to one valid fits-file
        extno          - int OR None  - Extension number of image data
                                      - (if multi-extension)
        verbose        - bool         - Print extra stuff if True
        EXPTIME        - str OR float - Exposure time in seconds
                                      - (keyword or value)
        MAGZERO        - str OR float - Zeropoint magnitude (keyword or value)
        MJD_START      - str OR float - MJD at start of exposure
                                      - (keyword or value)
        GAIN           - str OR float - Gain value (keyword or value)
        FILTER         - str          - Filter name (keyword)
        NAXIS1         - str OR int   - Number of pixels along axis 1
        NAXIS2         - str OR int   - Number of pixels along axis 2
        INSTRUMENT     - str          - Instrument name (keyword)

        A float/integer value can be defined for most keywords,
        rather than a keyword name; this will use that value
        rather than searching for the keyword in the headers.
        INSTRUMENT and FILTER obviously can't be floats/integers,
        so use a leading '-' to specify a value
        rather than a keyword name to use.

        output:
        -------
        self.header_keywords - dictionary - A bunch of header keyword names.
                                          - Needed later at all???
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
        self.header_keywords = {'EXPTIME': 'EXPTIME',      # Exposure time [s]
                                'MAGZERO': 'MAGZERO',      # Zeropoint mag
                                'MJD_START': 'MJD-OBS',    # MJD at start
                                'GAIN': 'GAIN',            # Gain value
                                'FILTER': 'FILTER',        # Filter name
                                'NAXIS1': 'NAXIS1',        # Pixels along axis1
                                'NAXIS2': 'NAXIS2',        # Pixels along axis2
                                'INSTRUMENT': 'INSTRUME',  # Instrument name
                                }
        # This is not a comprehensive list, just the ones that I expect
        # we'll need, particularly those that I know have different names
        # from different telescopes.

        # Do a loop over the kwargs and see if any header keywords need updating
        # (becasue they were supplied)
        for key, value_supplied in kwargs.items():
            if key in self.header_keywords:
                self.header_keywords[key] = value_supplied

        # Read the file:
        if extno is None:
            print('Warning: Treating this as a single extension file.')
            extno = 0

        with fits.open(filename) as han:
            self.data = han[extno].data
            self.header = han[extno].header  # Header for the extension
            self.header0 = han[0].header  # Overall header for mosaic, ext0

        # Read the WCS
        self.WCS = wcs.WCS(self.header)

        # Use the defined keywords to save the values into self.
        # Search both headers if neccessary (using keyValue)
        for key, use in self.header_keywords.items():
            if key not in ['INSTRUMENT', 'FILTER', 'NAXIS1', 'NAXIS2']:
                # Most keywords can just have a float value defined
                # instead of keyword name, that's what the if type is about.
                setattr(self, key, (float(self._find_key_value(use))
                                    if type(use) != float else use))
            elif key in ['NAXIS1', 'NAXIS2']:
                # Some keywords can have an integer value defined
                # instead of keyword name
                setattr(self, key, (int(self._find_key_value(use))
                                    if type(use) != int else use))
            elif key == 'INSTRUMENT':
                # INSTRUMENT and FILTER obviously can't be floats,
                # so use a leading '-' to specify a value to use
                # instead of a keyword name.
                self.INSTRUMENT = (self._find_key_value(use) if use[0] != '-'
                                   else use[1:])
            elif key == 'FILTER':
                # Filter only wants the first character of the supplied,
                # not all the junk (telescopes usually put numbers after)
                self.FILTER = (self._find_key_value(use)[0] if use[0] != '-'
                               else use[1])

        # Also define the middle of the exposure:
        self.MJD_MID = self.MJD_START + self.EXPTIME / 172800.0

        print('{}\n'.format((self.EXPTIME, self.MAGZERO, self.MJD_START,
                             self.MJD_MID, self.GAIN, self.FILTER,
                             self.NAXIS1, self.NAXIS2, self.INSTRUMENT))
              if verbose else '', end='')

    def _find_key_value(self, key):
        """
        First checks extension header for keyword; if fails, checks main header.
        This is neccessary because some telescopes put things like the EXPTIME
        in the main header, while others put it in the header for each
        extension and NOT the main one.

        input:
        ------
        key - str - keyword to look for

        output:
        -------
        value of keyword found in headers
        """
        try:
            value = self.header[key]
        except KeyError:
            value = self.header0[key]
        return value


# END
