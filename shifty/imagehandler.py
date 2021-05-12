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
from datetime import datetime
import copy
import numpy as np

from astropy.io import fits
from astropy import wcs
from astropy.nddata import CCDData
from ccdproc import wcs_project  # , Combiner


# -----------------------------------------------------------------------------
# Any local imports
# -----------------------------------------------------------------------------
sys.path.append(os.path.dirname(os.path.dirname(
                os.path.realpath(__file__))))
from shifty.downloader import Downloader
from shifty.known import Known


# -----------------------------------------------------------------------------
# Various class definitions for *data import * in shifty
# -----------------------------------------------------------------------------

class OneImage(Downloader):
    '''
    (1)Loads a fits-file,
       including data array, WCS, header and select header keywords
    (2)Calculate pixel coordinates? (Seems unneccessary at this stage?)

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


class DataEnsemble(Downloader):
    '''
    (1)Loads a list of fits-files
    (2)Stores image data in a cube array, WCS in an array, headers in an array.
    (3)Can reproject the data to a common WCS if neccesarry.

    methods:
    --------
    reproject_data()

    main public method:
    -------------------
    reproject_data()


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
        # - Allow ourselves to use Downloader methods
        Downloader.__init__(self,)
        self.filename = filename
        if filename is not None:
            # Set some values
            self.extno = extno
            self.reprojected = False
            self.read_fits_files(verbose=verbose, **kwargs)
        else:
            self.extno = None
            self.reprojected = None
            self.data = None
            self.WCS = None
            self.header = None

    def read_fits_files(self, verbose=False, **kwargs):
        '''
        Read in the files defined in self.filename
        from extension[s] in self.extno.

        input:
        ------
        self
        verbose        - bool         - Print extra stuff if True
        kwargs passed along into OneImage
        '''
        # Do some awesome stuff!!!
        datacube = []
        wcscube = []
        headercube = []
        for i, filei in enumerate(self.filename):
            exti = self.extno if isinstance(self.extno, int) else self.extno[i]
            OneIm = OneImage(filei, exti, verbose, **kwargs)
            datacube.append(OneIm.data)
            wcscube.append(OneIm.WCS)
            headercube.append(OneIm.header)
        self.data = np.array(datacube)
        self.WCS = np.array(wcscube)
        self.header = headercube

    def reproject_data(self, target=0):
        '''
        Reprojects each layer of self.data to be aligned, using their WCS.
        By default aligns everything else with the first image,
        but this can be changed by setting target.

        inputs:
        -------
        self.data, self.WCS, self.header
        target    - int - Index of the target image to align relative to

        outputs:
        --------
        self.data
        self.wcs
        self.header
        self.reprojected = True
        '''
        if self.reprojected:
            print('Data has already been aligned and reprojected! '
                  'Doing nothing!')
        else:
            for i, dat in enumerate(self.data):
                print(f"Reprojecting image {i}")
                if i != target:  # Don't reproject the target image, duh!
                    # Do the reprojection
                    reprojecti = wcs_project(CCDData(dat, wcs=self.WCS[i],
                                             unit='adu'), self.WCS[target])
                    # Append to list
                    self.data[i, :, :] = reprojecti.data
                    # Update WCS, both in self.wcs and self.header
                    self.WCS[i] = self.WCS[target]
                    del self.header[i]['CD?_?']  # required for update to work
                    self.header[i].update(self.WCS[i].to_header(relax=True))
                    # Add a comment to the header about the reprojection
                    now = str(datetime.today())[:19]
                    self.header[i]['COMMENT'] = (f'Data was reprojected to WCS '
                                                 f'of file {target} at {now}')
            self.reprojected = True


class DataHandler():
    '''
    (1) Loads a list of fits-files, via DataEnsemble
    (2) Shift data arrays by integer shifts
    (3) Stack the data in various ways
    (4) Save stacks
    (5) Interact with 'known'?
    (6) Interact with SourceExtractor?

    methods:
    --------
    integer_shift
    stack
    save_stack

    main public method:
    -------------------
    integer_shift
    stack

    '''
    # Turning off some stupid syntax-checker warnings:
    # pylint: disable=too-many-instance-attributes
    # Why on earth should an object only have 7 attributes?
    # That seems dumb. Turning off this warning.
    # pylint: disable=too-few-public-methods
    # I get why an object should have at least two public methods in order to
    # not be pointless, but it's an annoying warning during dev. Turning off.
    # pylint: disable=attribute-defined-outside-init
    # Why should a method not be allowed to redefine an attribute?
    # It doesn't make any sense.

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
        # Set some values
        self.filename = filename
        self.extno = extno
        if filename is not None:
            self.image_data = DataEnsemble(filename, extno, verbose, **kwargs)
        else:
            self.image_data = DataEnsemble()
        self.shifted_data = DataEnsemble()
        self.stacked_data = DataEnsemble()
        # Do some awesome stuff!!!

    def integer_shift(self, shifts, padmean=False):
        '''
        Shift some data.

        inputs:
        -------
        self.image_data
        shifts    - list of int OR array of int - Shape (N_images, 2)

        outputs:
        --------
        self.shifted_data - data array, updated WCS and header
        '''
        # Make sure shifts is an array of integers:
        shifts = np.array(shifts).round(0).astype(int)
        # Make sure shifts are all positive (by subtracting smallest value)
        shifts[:, 0] -= shifts[:, 0].min()
        shifts[:, 1] -= shifts[:, 1].min()
        # Find max shifts, to know how much to pad
        ymax = shifts[:, 0].max()
        xmax = shifts[:, 1].max()
        # Populate self.shifted_data.header and .WCS by copying
        # from self.image_data, then updating
        self.shifted_data.header = copy.deepcopy(self.image_data.header)
        self.shifted_data.WCS = copy.deepcopy(self.image_data.WCS)
        # Do the shifting in a loop over each image
        shifted = []
        for i, dat in enumerate(self.image_data.data):
            print(f'Shifting image {i} by {shifts[i]}')
            pad_value = np.mean(dat) if padmean else np.nan  # Mean or NaN
            pad_size = ((shifts[i, 0], ymax - shifts[i, 0]),  # Size of pad
                        (shifts[i, 1], xmax - shifts[i, 1]))  # on 4 sides
            # Shift the array by adding a padding around the edge
            shifted.append(np.pad(dat, pad_size, constant_values=pad_value))
            # Update WCS, both in .wcs and .header
            self.shifted_data.WCS[i].wcs.crpix += (shifts[i, 1], shifts[i, 0])
            del self.shifted_data.header[i]['CD?_?']  # required to update
            updated_WCS_keys = self.shifted_data.WCS[i].to_header(relax=True)
            self.shifted_data.header[i].update(updated_WCS_keys)
            # Add a comment to the header about the reprojection
            now = str(datetime.today())[:19]
            com_str = (f'At {now} data was padded by {pad_size}')
            self.shifted_data.header[i]['COMMENT'] = com_str
            com_str = (f'CRPIX keywords of WCS updated to reflect padding')
            self.shifted_data.header[i]['COMMENT'] = com_str
        # Turn the shifted arrays into a 3D array and stick into self.
        self.shifted_data.data = np.array(shifted)

    def shift_stack_by_name(self, object_name, obs_code,
                            object_type='smallbody',
                            padmean=False, **stack_args):
        '''
        Takes a known object (by name), and shifts the images according to the
        JPL Horizons ephemeris of the object.

        input:
        ------
        self.image_data
        object_name    - string - Name of object to be looked up in Horizons
        obs_code       - string - Observatory code
        padmean        - bool   - Whether to pad with the mean (True) or
                                  np.NaN (False)

        'obs_code' must be a Horizons obs code (eg. '500@-95' for Tess).
        Other optional keywords that selt.stack takes (such as 'which_WCS',
        'padmean', 'median_combine' and 'save_to_filename') are also accepted
        and passed along.

        output:
        -------
        self.shifted_data - data array, updated WCS and header
        '''
        shifts = self._calculate_shifts_from_known(object_name=object_name,
                                                   obs_code=obs_code,
                                                   object_type=object_type)
        self.integer_shift(shifts, padmean=padmean)
        print(shifts)
        self.stack(shifted=True, **stack_args)

    def stack(self, shifted=False, median_combine=False,
              save_to_filename='', which_WCS='middle'):
        '''
        Stacks the data.

        inputs:
        -------
        self.data or self.shifted_data
        shifted           - boolean        - Whether to use the shifted
                                           - or plain data
        median_combine    - boolean        - Use median if True, mean otherwise
        save_to_filename  - string or None - if string, saves to that file
        which_WCS         - string or int  - 'first', 'middle', 'last' or int
                                           - if int, it's used as the index

        outputs:
        --------
        self.stack

        A median stack is often a little deeper than the mean stack,
        but it does not preserve photon count in a photometric way.
        So median stack is good for finding objects, mean best for photometry.

        which_WCS defaults to using the WCS of the middle image.
        It does not actually matter, though.
        '''
        print('Combining images', end='')
        data = self.shifted_data.data if shifted else self.image_data.data
        sWCS = self.shifted_data.WCS if shifted else self.image_data.WCS
        # This method is much slower (factor 10-40) for some reason:
#        combiner = Combiner([CCDData(dati, unit='adu') for dati in data])
        if median_combine:  # slower; only if median is desired.
            print(' using median stacking.')
#            self.stacked_data.data = combiner.median_combine()  # Slow method
            self.stacked_data.data = np.median(data, 0)
        else:  # Default
            print(' using mean stacking')
#            self.stacked_data = combiner.average_combine()  # Slow method
            self.stacked_data.data = np.mean(data, 0)
        if isinstance(which_WCS, int):
            wcsidx = which_WCS
        elif isinstance(which_WCS, str):
            wcsidx = (0 if which_WCS.lower() == 'first'
                      else -1 if which_WCS.lower() == 'last'
                      else int(len(data) / 2))
        self.stacked_data.WCS = sWCS[wcsidx]
        if save_to_filename != '':
            self.save_stack(filename=save_to_filename)

    def save_stack(self, filename='stack.fits'):
        '''Save a stack to a fits file.'''
        save_fits(self.stacked_data, filename)

    def save_shifted(self, filename='shift'):
        '''Save the shifted images to fits files.'''
        save_fits(self.shifted_data, filename)

    def _calculate_shifts_from_known(self, **KnownArgs):
        '''
        Takes a known object and calculates its RA and Dec at the times
        in the FITs headers, converts to pixels and calculates the shift.

        input:
        ------
        Must have:
        obs_code    - string - Observatory code

        At least one of:
        object_name - string - Name of object to be looked up in Horizons
        orbit       - object - An orbit (not yet implemented)

        Either 'object_name' or 'orbit' must be supplied (if both are given,
        'orbit' is ignored). If 'object_name' is given, 'obs_code' must be a
        Horizons observatory code (eg. '500@-95' for Tess). If 'orbit' is given,
        'obs_code' must be ??? (not implemented yet, but probably MPC format).
        Other optional keywords that known.Known takes (such as 'object_type')
        are also accepted and passed along.

        output:
        -------
        shifts - Nx2 array - The calculated pixel shifts
        '''
        times = np.array([dh['SHIFTY_MJD_MID'] + 2400000.5 for dh
                          in self.image_data.header])
        if not (('object_name' in KnownArgs) or ('orbit' in KnownArgs)):
            raise TypeError('Keyword arguments "object_name" or "orbit" must '
                            'be supplied.')

        K_obj = Known(times=times, **KnownArgs)
        shifts = np.zeros([len(times), 2])
        try:  # try using all transformation parameters
            pix0 = self.image_data.WCS[0].all_world2pix(K_obj.RA[0],
                                                        K_obj.Dec[0], 0,
                                                        ra_dec_order=True)
            for i in np.arange(len(times)):
                pixi = self.image_data.WCS[i].all_world2pix(K_obj.RA[i],
                                                            K_obj.Dec[i], 0,
                                                            ra_dec_order=True)
                shifts[i, :] = np.flip(pix0) - np.flip(pixi)
        except wcs.NoConvergence:  # if fail, do basic WCS transform
            pix0 = self.image_data.WCS[0].wcs_world2pix(K_obj.RA[0],
                                                        K_obj.Dec[0], 0,
                                                        ra_dec_order=True)
            for i in np.arange(len(times)):
                pixi = self.image_data.WCS[i].wcs_world2pix(K_obj.RA[i],
                                                            K_obj.Dec[i], 0,
                                                            ra_dec_order=True)
                shifts[i, :] = np.flip(pix0) - np.flip(pixi)

        return shifts


# -------------------------------------------------------------------------
# These functions really don't need to be methods, and therefore aren't.
# No need to over-complicate things.
# -------------------------------------------------------------------------

def save_fits(DataEnsembleObject, filename='data.fits'):
    '''
    Save some data to [a] fits file[s].

    input:
    ------
    DataEnsembleObject - DataEnsemble Object - contains data and header.
    filename           - str                 - name or shared start of filename

    output:
    -------
    one or more fits file.
    If multiple, they are numbered *_000.fits, *_001.fits, etc.
    '''
    if len(DataEnsembleObject.data.shape) == 2:
        filename = filename if '.fits' in filename else filename + '.fits'
        hdu = fits.PrimaryHDU(data=DataEnsembleObject.data,
                              header=DataEnsembleObject.header)
        print(f'Saving to file {filename}')
        hdu.writeto(filename, overwrite=True)
    elif len(DataEnsembleObject.data.shape) == 3:
        for i, data in enumerate(DataEnsembleObject.data):
            filenamei = filename.replace('.fits', '')
            filenamei += f'_{i:03.0f}.fits'
            hdu = fits.PrimaryHDU(data=data,
                                  header=DataEnsembleObject.header[i])
            print(f'Saving to file {filenamei}')
            hdu.writeto(filenamei, overwrite=True)
    else:
        raise ValueError('The input is not a valid DataEnsemble object.')

    print('Done!')


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
    header_comments = {'EXPTIME': 'Exposure time [s]',
                       'MAGZERO': 'Zeropoint magnitude',
                       'MJD_START': 'MJD at start of exposure',
                       'GAIN': 'Gain value',
                       'FILTER': 'Filter letter',
                       'NAXIS1': 'Pixels along axis1',
                       'NAXIS2': 'Pixels along axis2',
                       'INSTRUMENT': 'Instrument name',
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
        header[f'SHIFTY_{key}'] = (key_values[key], header_comments[key])
        header['COMMENT'] = (f'SHIFTY_{key} added by SHIFTY'
                             f' at {str(datetime.today())[:19]}')
        header['COMMENT'] = f'SHIFTY_{key} derived from {use}'

    # Also define the middle of the exposure:
    key_values['MJD_MID'] = (key_values['MJD_START'] +
                             key_values['EXPTIME'] / 172800.0)
    header[f'SHIFTY_MJD_MID'] = (key_values['MJD_MID'],
                                 'MJD at middle of exposure')
    header['COMMENT'] = (f'SHIFTY_MJD_MID added by SHIFTY'
                         f' at {str(datetime.today())[:19]}')
    header['COMMENT'] = (f'SHIFTY_MJD_MID derived from '
                         f'{header_keywords["MJD_START"]} '
                         f'and {header_keywords["EXPTIME"]}')

    print('{}\n'.format((key_values)) if verbose else '', end='')
    return data, header, header0, wcs.WCS(header), header_keywords, key_values


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
