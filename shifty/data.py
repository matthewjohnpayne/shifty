

# -------------------------------------------------------------------------------------
# Third party imports
# -------------------------------------------------------------------------------------
import os, sys
import astropy
from astropy.io import fits
from collections import OrderedDict
import numpy as np

# -------------------------------------------------------------------------------------
# Any local imports
# -------------------------------------------------------------------------------------
# N/A


# -------------------------------------------------------------------------------------
# Various class definitions for *Data Handling* in shifty
# -------------------------------------------------------------------------------------

class ImageDataSet():
    '''
        A set of images w/ times & wcs which are suitable for stacking (e.g. stars have been masked, bad cadences removed, …)
        
        Inputs:
        -------
        (i) images:
        list of `astropy.io.fits.hdu.hdulist.HDUList` objects
        
        (ii) obs_code:
        code (string) to uniquely specify observatory, and hence allow determination of observatory position as a func of time
        
        Methods:
        --------
        observatory_position (JPL code)
        get_observatory_barycentric_positions()
        get_theta_wcs()  expose the transformation of pixel to theta space

    '''

    def __init__(self, HDUs , obs_code) :

        # Each HDU as read from fits file will have several components
        # - In the TESS files there are 3xHeaders (Primary, Calibrated-Data, Uncertainty)
        # - In the TESS files there are 2xData (Calibrated-Data, Uncertainty)
        # Do we want to split them out at this point ?
        self.images     = images
        self.obs_code   = obs_code


    def generate_observatory_barycentric_positions(self):
        '''list_of_barycentric_positions_one_for_each_image'''
        return [calculate_barycentric_position(image.header['obs_date'] , self.obs_code) for image in self.images]

    def generate_theta_coordinates(self):
        ''' 
            Each pixel in an image corresponds to some RA,Dec coord (can get from wcs)
             - N.B. an (RA,Dec) can be represented as a unit-vector, u=(u_x, u_y, u_z)
             
            We want to rotate/transform the coordinates of each pixel in each image into a single, standard, set of projection-plane coordinates
        
            We are *NOT* touching flux values
            
            For each image this will return a 2D array of "theta-coordinates", one for each pixel in the image
            
            We need to think about whether/how the tangent plane is pixelized

        '''
        pass



class ImageLoader():
    '''
        Parent class for ...
         - TessImageLoader, HubbleImageLoader, PanstarrsImageLoader, ...
        
        * inputs:
        sector number, camera, chip, detector_range, bad_data_thresholds
        
        * methods
        _load_images()
        _remove_stars()  (i.e. set pixels that contain stars to NaN)
        _remove_bad_cadences()
        _remove_scattered_light_problem_areas()
        _remove_strap_regions()
        
        * main public method:
        get_image_data_set() => Returns an ImageDataSet
        

    '''
    
    def __init__(self, ) :
        
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
    # Data directories / Storage / etc
    # -------------------------------------------------------------------------------------
    def _fetch_data_directory(self):
        '''
            Returns the default path to the directory where data will be downloaded.
        
            By default, this method will return ~/.shifty/data
            and create this directory if it does not exist.
        
            If the directory cannot be accessed or created, then it returns the local directory (".")
            
            Returns
            -------
            data_dir : str
                Path to location of `data_dir` where data (FITs files) will be downloaded
        '''

        data_dir = os.path.join(os.path.expanduser('~'), '.shifty_data')
        
        # If it doesn't exist, make a new data directory
        if not os.path.isdir(data_dir):
            
            try:
                os.mkdir(data_dir)
            
            # downloads locally if OS error occurs
            except OSError:
                warnings.warn('Warning: unable to create {}. '
                              'Download directory set to be the current working directory instead.'.format(data_dir))
                data_dir = '.'
                                                
        return data_dir
    

    # -------------------------------------------------------------------------------------
    # The methods below are for the loading of *general* fits-format data files
    # -------------------------------------------------------------------------------------

    def _load_images(self, file_spec_container ):
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
        
        # if file_spec_container happens to be a single valid fits file(s), then we are good to go ...
        if isinstance(file_spec_container, str) and os.path.isfile(file_spec_container) :
            fits_file_paths = [file_spec_container]
        
        # if file_spec_container happens to be an iterable containing only valid fits file(s), then we are good to go ...
        elif isinstance(file_spec_container, (list, tuple, dict) ) and np.all( [ os.path.isfile(_) for _ in file_spec_container] ):
            fits_file_paths = [ _ for _ in file_spec_container]
        
        # no other handling-methods currently in place (but CHILD can attempt some PRE-FILTERING)
        else:
            print('file_spec_container in *_load_images()* is of type [%r] and has the content ... \n\t %r ' % ( type(file_spec_container), file_spec_container) )
            print(' *** At present, no method is in place to interpret this input *** ')
            print(' ******                No files will be loaded               ******')
            fits_file_paths = []
        
        # open filepaths
        return [ self._load_image(fp) for fp in fits_file_paths ]
    
    def _load_image(self , fits_file_path):
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
        return fits.open(fits_file_path) if os.path.isfile(fits_file_path) else None








class TESSImageLoader(ImageLoader):
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
        
        # - Allow ourselves to use ImageLoader methods
        ImageLoader.__init__(self, )
    
        # - Define some important TESS-related quantities
        self.file_spec_keywords = ['sector', 'camera', 'chip']
        self.testing_keywords   = ['DEV', 'TEST']
        self.obs_code            = 'C57'
    
    # -------------------------------------------------------------------------------------
    # Public Methods
    # -------------------------------------------------------------------------------------
    def get_image_data_set(self, file_spec_container={}, cleaning_parameter_dict ={} ):
        '''
            Overall image loader as envisaged by Geert & Matt over coffee
            - Gets data from file(s)
            - Does "cleaning"
            - Creates ImageDataSet object
            
            Input:
            ------
            (1) cleaning_spec_dict
             - container with enough info to specify what files to open
             - is passed-through-to and interpreted-by *_load_images()*
             - [[probably a cleaner way to do this involving **kwargs, but I just want to make progress]]
            (2) cleaning_spec_dict
             - option args to pass through to various sub-functions
            
            Returns:
            --------
            ImageDataSet
            '''
        
        # Load the data from files
        HDUs = self._load_images( file_spec_container )
        
        # clean the data
        # - [[could just ovewrite original HDUs with clean version: will leave both for now]]
        cleanHDUs = self._clean_data(HDUs, cleaning_parameter_dict )
        
        # create & return ImageDataSet
        return ImageDataSet(cleanHDUs, self.obs_code)

    # -------------------------------------------------------------------------------------
    # The methods below are for the "CLEANING" of TESS data
    # -------------------------------------------------------------------------------------
    def _clean_data(self, HDUs, cleaning_parameters):
        ''' 
            Wrapper around the various "cleaning" methods below
            - Simply provided as a means to enable simple toggling on/off of functions 
            - Uses **kwargs to control what is/isn't evaluated
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
        # [[ NOTICE THE IMPLICIT DESIGN CHOICE THAT ALL CLEANING FUNCTIONS MUST RETURN HDUs]]
        for key, func_to_run in cleaning_function_dict.items():
            # run a function if it is included as True in cleaning_parameters (e.g. {'mask':True})
            if key in cleaning_parameters and cleaning_parameters[key]:
                HDUs = func_to_run(HDUs, cleaning_parameters)


        

    def _mask_stars(self,HDUs, **kwargs):
        ''' 
            We want to remove stars in some way
            Barentsen & Payne discussed a simple mask: i.e. set pixels that contain stars to NaN
            This would be done based on GAIA positions
            
            This is *NOT* subtraction (see _subtract_stars below )
            
            Presumably only one of _mask_stars / _subtract_stars is required, but I am 100% certain that Holman will at least want to experiment with subtraction
            
        '''
        return HDUs
    
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
    # The methods below relate to the aquisition of TESS data
    # -------------------------------------------------------------------------------------
    def _ensure_tess_subdirectory_structure_exists(self, sectorNumber, cameraNumber, chipNumber):
        '''
            Ensure subdirectory structure exists
             - Subdirectories organized as self.local_dir/sector/camera/chip
            '''
        try:
            assert os.path.isdir(self.local_dir)
            
            subDir = os.path.join(self.local_dir, sectorNumber )
            if not os.path.isdir(subDir): os.mkdir(subDir)
            
            subDir = os.path.join(subDir, cameraNumber )
            if not os.path.isdir(subDir): os.mkdir(subDir)
            
            subDir = os.path.join(subDir, chipNumber )
            if not os.path.isdir(subDir): os.mkdir(subDir)
            
            assert os.path.isdir(subDir)
            return True
        except OSError:
            sys.exit('Problem with sub-directory creation in ..._ensure_tess_subdirectory_structure_exists()...')
    

    def _fetch_tess_fits_filepath(self, fit_filename):
        '''
            Parses the name of a tess fits-file and decides on a filepath
             - fits filename looks like tess2018292175940-s0004-1-2-0124-s_ffic.fits
            Filepath is just subdirectories within self.local_dir
            Subdirectories organized as self.local_dir/sector/camera/chip
        '''
        # can split on "-" to extract sector number, camera & chip
        sectorNumber   = int(fit_filename.split("-")[1][-4:])
        cameraNumber   = int(fit_filename.split("-")[2])
        chipNumber     = int(fit_filename.split("-")[3])
    
        # derive filepath
        tessDir     = os.path.join(self.local_dir, str(sectorNumber), str(cameraNumber), str(chipNumber) )
        filepath    = os.path.join(tessDir, fit_filename)
        
        # ensure directory-structure exists
        self._ensure_tess_subdirectory_structure_exists(str(sectorNumber), str(cameraNumber), str(chipNumber))
        
        return {
            'sectorNumber'  : sectorNumber,
            'cameraNumber'  : cameraNumber,
            'chipNumber'    : chipNumber,
            'filepath'      : filepath
        }

    def _download_download_script(self , sectorNumber ):
        '''
            There are scripts available online to facilitate bulk downloads
            E.g. https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_4_ffic.sh
            
            This function will get the script that is relevant for a specific TESS sector
        '''
        filename    = 'tesscurl_sector_%s_ffic.sh' % sectorNumber
        outFilepath = os.path.join(self.local_dir , filename)
        command     = "curl -C - -L -o %s https://archive.stsci.edu/missions/tess/download_scripts/sector/%s" % ( outFilepath , filename)
        os.system(command)
        return outFilepath
        
    def _parse_download_script(self , shFile ):
        '''
            Parse a script downloaded using '_download_download_script'
             - I.e. read it to understand the names of the files that will be downloaed and hence where we want to file them away
            
            Contents look like 16,000 rows of ...
            curl -C - -L -o tess2018292175940-s0004-1-2-0124-s_ffic.fits https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/tess2018292175940-s0004-1-2-0124-s_ffic.fits
            
            Will output a dictionary containing ...
             - list of fits files
            E.g. [ 'tess2018292095940-s0004-1-1-0124-s_ffic.fits' , '']
            
             - Will also parse the fits files to understand their sector, camera, chip, ...
            E.g. [(4,1,1,)]
            
             - Will also regenerate a 'corrected' curlcommand to put the file into the location we want 
            E.g. 'curl -C - -L -o /Users/matthewjohnpayne/.shifty_data/4/4/4/tess2018292105940-s0004-4-4-0124-s_ffic.fits https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/tess2018292105940-s0004-4-4-0124-s_ffic.fits'
        '''
        # read contents
        with open(shFile, 'r') as fh:
            data = [_ for _ in fh.readlines() if len(_.split()) > 5 and _[0]!='#']
        
        # fits files are in the fifth secn if we split on white-space
        returnDict = {  'fits_files':[_.split()[5] for _ in data],
                        'sectorNumbers':[], 'cameraNumbers':[], 'chipNumbers':[], 'filepaths':[], 'curlCommands':[] }
        
        # can then use _fetch_tess_fits_filepath() to get the sector number, camera, chip & destination-filepath
        for ff, originalCurl in zip(returnDict['fits_files'], data):
            result = self._fetch_tess_fits_filepath(ff)
        
            # just getting data out of the returned dictionary
            returnDict['sectorNumbers'].append(result['sectorNumber'])
            returnDict['cameraNumbers'].append(result['cameraNumber'])
            returnDict['chipNumbers'].append(result['chipNumber'])
            returnDict['filepaths'].append(result['filepath'])
        
            # edit the curl command to redirect output to the sub-directory that we want  ...
            returnDict['curlCommands'].append(' '.join( (*originalCurl.split()[:5], result['filepath'], *originalCurl.split()[6:] ) ) )
        
        return returnDict



    # -------------------------------------------------------------------------------------
    # The method(s) below are for the loading of TESS fits-format data files
    # ----------------------------------------------------- --------------------------------

    def _load_images(self, file_spec_container ):
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
        # If 'dev'/'test' is specified, then get limited test/development data-set
        # - Will match either STRING or DICTIONARY-ENTRIES
        # - N.B. self.testing_keywords ~ ['DEV', 'TEST']
        if  isinstance (file_spec_container, str)  \
            and np.any( [ True for _ in self.testing_keywords if _ in file_spec_container.upper() ] ) \
            or  \
            isinstance (file_spec_container, dict) \
            and np.any( [ True for _ in self.testing_keywords if _ in file_spec_container and file_spec_container[_]] ):

            # Use the test-data method to ensure files available locally
            fits_file_paths = self._ensure_test_data_available_locally()
            to_parent       = fits_file_paths
            
            

        # If the file_spec_container contains the keywords required to specify ...
        # ... a set of TESS fits-files, then get the required filepaths
        # - N.B. self.file_spec_keywords = ['sector', 'camera', 'chip']
        #
        # [[ *** AS WRITTEN THIS WILL ONLY RETURN THE DATA FOR A SINGLE CHIP *** ]]
        #
        elif    isinstance (file_spec_container, dict) \
            and np.all( [ True in file_spec_container for _ in self.file_spec_keywords] ):
                directory_path = os.path.join( self.local_dir, *[ '%s' % str(file_spec_container[_]) for _ in self.file_spec_keywords] )
                to_parent      = glob.glob( os.path.join(directory_path , '*.fits'))

        # Might be something that the parent knows how to parse ...
        else:
            to_parent      = file_spec_container
                
        # Hand off to PARENT to do the rest of the loading
        return super(TESSImageLoader, self)._load_images( to_parent )


    # -------------------------------------------------------------------------------------
    # The methods below are for the loading of a limited, pre-defined set of TESS data
    # Intended for testing & development
    # -------------------------------------------------------------------------------------

    """
    def _load_test_data(self,):
        '''
            Load limited set of fits files for testing from local disk
            
            As per the overall expectation of ImageLoader/TESSImageLoader/get_image_data_set, this needs to return an "ImageDataSet"
        '''
        # If the test data doesn't exist on local disk, go and download the data from MAST
        self._ensure_test_data_available_locally()
    
        # Now that we are guaranteed that the data exists, we can open it
        # *** MJP: I Expect that the subsequent command/call will be updated (at some time in the future) to use a more general load function
        return self._load_images( self._define_test_data() )
    """
    

    def _ensure_test_data_available_locally(self,):
        '''
            Check whether local versions of the fits files are available on local disk
            If they do not exist, download from MAST
        '''
        for fits_file in self._define_test_data():
            if not os.path.isfile(  self._fetch_tess_fits_filepath( fits_file )['filepath'] ):
                self._download_test_data(fits_file)
        return self._define_test_data()

    def _download_test_data(self , fits_file ):
        '''
            Download fits file from MAST to local disk
            Note that file is saved to specific directory, "local_dir"
        '''
        command = "curl -C - -L -o %s https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/%s" % ( os.path.join(self.local_dir, fits_file) , fits_file)
        os.system(command)
        assert os.path.isfile( os.path.join(self.local_dir, fits_file) )


    def _define_test_data(self,):
        '''
           Define the fits files which we want to use for testing
           These were selected as being the first 10 files from sector-4 (see tesscurl_sector_4_ffic.sh on ...)
           Sector 4 was selected because many early teething-troubles (with pointing, etc) were overcome by this time
           
        '''
        return [ self._fetch_tess_fits_filepath(_)['filepath'] for _ in \
                    [
                    "tess2018292095940-s0004-1-1-0124-s_ffic.fits",
                    "tess2018292102940-s0004-1-1-0124-s_ffic.fits",
                    "tess2018292105940-s0004-1-1-0124-s_ffic.fits",
                    "tess2018292112940-s0004-1-1-0124-s_ffic.fits",
                    "tess2018292115940-s0004-1-1-0124-s_ffic.fits",
                    "tess2018292122940-s0004-1-1-0124-s_ffic.fits",
                    "tess2018292125940-s0004-1-1-0124-s_ffic.fits",
                    "tess2018292132940-s0004-1-1-0124-s_ffic.fits",
                    "tess2018292135940-s0004-1-1-0124-s_ffic.fits",
                    "tess2018292142940-s0004-1-1-0124-s_ffic.fits"
                    ]
                ]





class HSTImageLoader(ImageLoader):
    '''
        You know we'll want to do it at some point !!!
    '''
    pass

