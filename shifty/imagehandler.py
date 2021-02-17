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
    (1)Loads a fits-file
    (2)Reads WCS
    (3)Calculate pixel coordinates
    (4)Interact with 'known'?

    methods:
    --------
    loadImageAndHeader()
    others?

    main public method:
    -------------------
    loadImageAndHeader()


    '''
    # Turning off some stupid syntax-checker warnings:
    # pylint: disable=too-many-instance-attributes
    # Why on earth should an object only have 7 attributes?
    # That seems dumb. Turning off this warning.
    # pylint: disable=too-few-public-methods
    # I get why an object should have at least two public methods in order to
    # not be pointless, but it's an annoying warning during dev. Turning off.

    def __init__(self, filename=None, extno=0, verbose=False, **kwargs):
        '''
        inputs:
        -------
        filename       - str          - filepath to one valid fits-file
        extno          - int          - Extension number of image data
                                      - (0 if single-extension)
        verbose        - bool         - Print extra stuff if True
        EXPTIME        - str OR float - Exposure time in seconds
                                      - (keyword or value)
        MAGZERO        - str OR float - Zeropoint magnitude
                                      - (keyword or value)
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
        '''
        # - Allow ourselves to use Downloader methods
        Downloader.__init__(self,)
        # Make readOneImageAndHeader a method even though it doesn't need to be
        self.readOneImageAndHeader = readOneImageAndHeader
        # - Local directory for saving data
        self.local_dir = self._fetch_data_directory()
        # Initialize some attributes that will get filled later
        self.key_values = {}       # filled out by loadImageAndHeader
        self.header_keywords = {}  # filled out by loadImageAndHeader
        self.WCS = None            # filled out by loadImageAndHeader
        self.header = None         # filled out by loadImageAndHeader
        self.header0 = None        # filled out by loadImageAndHeader
        self.data = None           # filled out by loadImageAndHeader
        self.filename = filename
        self.extno = extno
        if filename:
            (self.data, self.header, self.header0, self.WCS,
             self.header_keywords, self.key_values
             ) = self.readOneImageAndHeader(filename, extno,
                                            verbose=verbose, **kwargs)

    # -------------------------------------------------------------------------
    # Public Methods
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # The methods below are for the loading of *general* fits-format data files
    # -------------------------------------------------------------------------


class ImageEnsemble(OneImage):
    '''
    (1)Loads a list of fits-files
    (2)Reads WCS
    (3)Calculate pixel coordinates
    (4)Interact with 'known'?

    methods:
    --------
    loadImageAndHeader()
    others?

    main public method:
    -------------------
    loadImageAndHeader()


    '''
    # Turning off some stupid syntax-checker warnings:
    # pylint: disable=too-few-public-methods
    # I get why an object should have at least two public methods in order to
    # not be pointless, but it's an annoying warning during dev. Turning off.

    def __init__(self, filename=None, extno=0, verbose=False, **kwargs):
        '''
        inputs:
        -------
        filename       - list of str  - list of filepaths to valid fits-files
        extno          - int OR       - Extension number to use for all images
                         list of int  - list of extension to use for each image
        verbose        - bool         - Print extra stuff if True
        EXPTIME        - str OR float - Exposure time in seconds
                                      - (keyword or value)
        MAGZERO        - str OR float - Zeropoint magnitude
                                      - (keyword or value)
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
        Not yet supported: list of a keyword values/names for each image, which
        hopefully should only ever be useful if attempting to stack images
        from different instruments; so low priority.
        '''
        # - Allow ourselves to use OneImage methods
        # This sets up readOneImageAndHeader as a method (even though it
        # doesn't have to be) and sets up the local_dir for saving data.
        OneImage.__init__(self,)
        # Set some values
        self.filename = filename
        self.extno = extno
        # Do some awesome stuff!!!
        datacube = []
        wcscube = []
        mjdcube = []
        for i, filei in enumerate(filename):
            exti = extno if isinstance(extno, int) else extno[i]
            OneIm = OneImage(filei, exti, verbose, **kwargs)
            datacube.append(OneIm.data)
            wcscube.append(OneIm.WCS)
            mjdcube.append(OneIm.key_values['MJD_MID'])
        self.data = np.array(datacube)
        self.WCS = np.array(wcscube)
        self.MJD = np.array(mjdcube)


# -------------------------------------------------------------------------
# These functions really don't need to be methods, and therefore aren't.
# No need to over-complicate things.
# -------------------------------------------------------------------------

def readOneImageAndHeader(filename=None, extno=0, verbose=False,
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
    data            - np.array   - the float array of pixel data
                                 - shape == (NAXIS2, NAXIS1)
    header          - fits.header.Header - Sort of like a dictionary
    header0         - fits.header.Header - Sort of like a dictionary
    WCS             - wcs.wcs.WCS - World Coordinate System plate solution
    header_keywords - dictionary - A bunch of header keyword names.
                                 - Needed later at all??? Not sure
    key_values      - dictionary - A bunch of important values

    Content of key_values:
    EXPTIME         - float      - Exposure time in seconds
    MAGZERO         - float      - Zeropoint magnitude
    MJD_START       - float      - MJD at start of exposure
    MJD_MID         - float      - MJD at centre of exposure
    GAIN            - float      - Gain value
    FILTER          - str        - Filter name
    NAXIS1          - int        - Number of pixels along axis 1
    NAXIS2          - int        - Number of pixels along axis 2
    '''

    # Check whether a filename is supplied.
    if filename is None:
        raise TypeError('filename must be supplied!')

    # Define default keyword names
    header_keywords = {'EXPTIME': 'EXPTIME',      # Exposure time [s]
                       'MAGZERO': 'MAGZERO',      # Zeropoint mag
                       'MJD_START': 'MJD-STR',    # MJD at start
                       'GAIN': 'GAINEFF',         # Gain value
                       'FILTER': 'FILTER',        # Filter name
                       'NAXIS1': 'NAXIS1',        # Pixels along axis1
                       'NAXIS2': 'NAXIS2',        # Pixels along axis2
                       'INSTRUMENT': 'INSTRUME',  # Instrument name
                       }
    key_values = {}

    # Do a loop over the kwargs and see if any header keywords need updating
    # (because they were supplied)
    for key, non_default_name in kwargs.items():
        if key in header_keywords:
            header_keywords[key] = non_default_name

    # Read the file. Do inside a "with ... as ..." to auto close file after
    with fits.open(filename) as han:
        data = han[extno].data
        header = han[extno].header  # Header for the extension
        # Overall header for whole mosaic, etx0:
        header0 = han[0].header  # pylint: disable=E1101 # Pylint stupid errors

    # Read the WCS
    WCS = wcs.WCS(header)

    # Use the defined keywords to save the values into key_values.
    # Search both headers if neccessary (using keyValue)
    for key, use in header_keywords.items():
        if key not in ['INSTRUMENT', 'FILTER', 'NAXIS1', 'NAXIS2']:
            # Most keywords can just have a float value defined
            # instead of keyword name, that's what the if type is about.
            key_values[key] = (use if isinstance(use, float) else
                               float(_find_key_value(header, header0, use)))
        elif key in ['NAXIS1', 'NAXIS2']:
            # Some keywords can have an integer value defined
            # instead of keyword name
            key_values[key] = (use if isinstance(use, int) else
                               int(_find_key_value(header, header0, use)))
        elif key == 'INSTRUMENT':
            # INSTRUMENT and FILTER obviously can't be floats,
            # so use a leading '-' to specify a value to use
            # instead of a keyword name.
            key_values[key] = (use[1:] if use[0] == '-'
                               else _find_key_value(header, header0, use))
        elif key == 'FILTER':
            # Filter only wants the first character of the supplied,
            # not all the junk (telescopes usually put numbers after)
            key_values[key] = (use[1] if use[0] == '-' else
                               _find_key_value(header, header0, use)[0])

    # Also define the middle of the exposure:
    key_values['MJD_MID'] = (key_values['MJD_START'] +
                             key_values['EXPTIME'] / 172800.0)

    print('{}\n'.format((key_values)) if verbose else '', end='')
    return data, header, header0, WCS, header_keywords, key_values


def _find_key_value(header1, header2, key):
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
        value = header1[key]
    except KeyError:
        value = header2[key]
    return value


# END
