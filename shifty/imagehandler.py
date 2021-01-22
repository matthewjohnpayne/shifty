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

class OneImage(Downloader):
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

    def __init__(self, filename=None, extno=0, verbose=False, **kwargs):
        # - Allow ourselves to use Downloader methods
        Downloader.__init__(self,)
        # - Local directory for saving data
        self.local_dir = self._fetch_data_directory()
        # Initialization of default standard Fits header keywords
        # - These may be overwritten by kwargs
        self.header_keywords = {'EXPTIME': 'EXPTIME',      # Exposure time [s]
                                'MAGZERO': 'MAGZERO',      # Zeropoint mag
                                'MJD_START': 'MJD_STR',    # MJD at start
                                'GAIN': 'GAINEFF',         # Gain value
                                'FILTER': 'FILTER',        # Filter name
                                'NAXIS1': 'NAXIS1',        # Pixels along axis1
                                'NAXIS2': 'NAXIS2',        # Pixels along axis2
                                'INSTRUMENT': 'INSTRUME',  # Instrument name
                                }
        # This is not a comprehensive list, just the ones that I expect
        # we'll need, particularly those that I know have different names
        # from different telescopes.
        self.key_values = {}  # filled out by loadImageAndHeader
        self.filename = filename
        self.extno = extno
        if filename:
            self.readOneImageAndHeader(filename, extno,
                                       verbose=verbose, **kwargs)

    # -------------------------------------------------------------------------
    # Public Methods
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # The methods below are for the loading of *general* fits-format data files
    # -------------------------------------------------------------------------

    def readOneImageAndHeader(self, filename=None, extno=0, verbose=False,
                              **kwargs):
        '''
        Reads in a fits file (or a given extension of one).
        Returns the image data, the header, and a few useful keyword values.

        input:
        ------
        filename       - str          - valid filepath to one valid fits-file
        extno          - int          - Extension number of image data
                                      - (0 if single-extension)
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
        self.key_values      - dictionary - A bunch of important values

        Content of self.key_values:
        EXPTIME         - float      - Exposure time in seconds
        MAGZERO         - float      - Zeropoint magnitude
        MJD_START       - float      - MJD at start of exposure
        MJD_MID         - float      - MJD at centre of exposure
        GAIN            - float      - Gain value
        FILTER          - str        - Filter name
        NAXIS1          - int        - Number of pixels along axis 1
        NAXIS2          - int        - Number of pixels along axis 2
        '''

        # If filename is supplied, use supplied filename and extension.
        # Otherwise, use those in self
        if (filename is None) & (self.filename is None):
            raise TypeError('filename must be supplied!')
        elif filename is not None:
            self.filename = filename
            self.extno = extno
        else:
            filename = self.filename
            extno = self.extno

        # Do a loop over the kwargs and see if any header keywords need updating
        # (becasue they were supplied)
        for key, non_default_name in kwargs.items():
            if key in self.header_keywords:
                self.header_keywords[key] = non_default_name

        # Read the file. Do inside a "with ... as ..." to auto close file after
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
                self.key_values[key] = (use if isinstance(use, float)
                                        else float(self._find_key_value(use)))
            elif key in ['NAXIS1', 'NAXIS2']:
                # Some keywords can have an integer value defined
                # instead of keyword name
                self.key_values[key] = (use if isinstance(use, int)
                                        else int(self._find_key_value(use)))
            elif key == 'INSTRUMENT':
                # INSTRUMENT and FILTER obviously can't be floats,
                # so use a leading '-' to specify a value to use
                # instead of a keyword name.
                self.key_values[key] = (use[1:] if use[0] == '-'
                                        else self._find_key_value(use))
            elif key == 'FILTER':
                # Filter only wants the first character of the supplied,
                # not all the junk (telescopes usually put numbers after)
                self.key_values[key] = (use[1] if use[0] == '-'
                                        else self._find_key_value(use)[0])

        # Also define the middle of the exposure:
        self.key_values['MJD_MID'] = (self.key_values['MJD_START'] +
                                      self.key_values['EXPTIME'] / 172800.0)

        print('{}\n'.format((self.key_values)) if verbose else '', end='')

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
