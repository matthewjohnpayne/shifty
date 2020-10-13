'''
   Classes / Methods for downloading fits files (e.g. from MAST) to the user's local disk
   
   Overengineed for something that won't be used much ...
'''

# -------------------------------------------------------------------------------------
# Third party imports
# -------------------------------------------------------------------------------------
import os, sys
import astropy
from collections import OrderedDict
import numpy as np
import functools
import glob

from astropy.io import fits
from astropy.time import Time

# -------------------------------------------------------------------------------------
# Any local imports
# -------------------------------------------------------------------------------------
# N/A


# -------------------------------------------------------------------------------------
# Various class definitions for *Data Downloading* in shifty
# -------------------------------------------------------------------------------------

class Downloader():
    '''
        Parent class for ...
         - TessDownloader, HSTDownloader, ...
         
        This parent class is NOT likely to be used directly by the user
         - The child classes (TessDownloader, HSTDownloader, ...) will be used instead

    '''
    
    def __init__(self, ) :
        
        # - Local directory for saving data
        self.local_dir = self._fetch_data_directory()

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
    # ...
    # -------------------------------------------------------------------------------------





class TESSDownloader(Downloader):
    '''
        Methods to download TESS data (from mast)
        Includes methods for downloading
        - FFIs
        - PRFs
        
    '''
        
    def __init__(self, ) :
        
        # - Allow ourselves to use any Downloader methods
        Downloader.__init__(self, )
        
        # - Local directory for saving data
        self.tess_dir = self._fetch_tess_data_directory()
    
    
    # -------------------------------------------------------------------------------------
    # Public methods for data download
    # -------------------------------------------------------------------------------------
    def download_sector(self, sectorNumber, cameraNumber, chipNumber):
        '''
            Download all fits files for a single sector
        '''
        for i in range(4):
            self.download_camera(sectorNumber, cameraNumber)
    
    def download_camera(self, sectorNumber, cameraNumber, chipNumber):
        '''
            Download all fits files for a single camera in a single sector
        '''
        for i in range(4):
            self.download_chip(sectorNumber, cameraNumber, chipNumber)

    def download_chip(self, sectorNumber, cameraNumber, chipNumber):
        '''
            Download all fits for a single chip in a single camera in a single sector
            Expect that thet will be put into directory like ... 
            ... ~/.shifty_data/tess/4/4/4/
        '''
        print('\t Downloading sectorNumber=%d, cameraNumber=%d, chipNumber=%d' % \
            (sectorNumber, cameraNumber, chipNumber))
            
        # Get a description of what files available for download
        returnDict = self._parse_download_script( self._get_download_script( sectorNumber ) )

        # Only execute the curl-commands appropriate for the sector/camera/chip
        for sec,cam,chp,fp, curl in zip(returnDict['sectorNumbers'],
                                        returnDict['cameraNumbers'],
                                        returnDict['chipNumbers'],
                                        returnDict['filepaths'],
                                        returnDict['curlCommands']):
            if sec == sectorNumber \
            and cam == cameraNumber \
            and chp == chipNumber \
            and not os.path.isfile(fp):
                try:
                    os.system(curl)
                except (Exception, OSError) as error :
                    print('sectorNumber, cameraNumber, chipNumber = ', sectorNumber, cameraNumber, chipNumber)
                    print('could not execute\n%r' % curl)
                    print(error)
                    print('will attempt to continue downloading other files ... ')
                    print()
                        
    def download_prf(self, ):
        '''
            download pixel response function models
             - from  https://archive.stsci.edu/missions/tess/models/prf_fitsfiles/
             
            for documentation & description:
             - see https://outerspace.stsci.edu/display/TESS/2.0+-+Data+Product+Overview#id-2.0-DataProductOverview-PixelResponseFunctions
        '''
        # At the end we want to pass back a list of prf filepaths
        prf_filepaths = []
        
        for sectorString in ['start_s0001','start_s0004']:
            
            # ensure local directory exists ...
            sector_directory = os.path.join( self._fetch_tess_prf_directory(), sectorString)
            if not os.path.isdir(sector_directory):
                os.mkdir(sector_directory)
            
            # get tar-files for sector
            for camera in range(1,5):
                
                # local filepath
                local_tar_filepath =  os.path.join( sector_directory, 'tess_prf_camera_%s.tgz' % str(camera) )
                
                # curl download command
                command = "curl -C - -L -o %s https://archive.stsci.edu/missions/tess/models/prf_fitsfiles/%s/tess_prf_camera_%s.tgz" % \
                    (local_tar_filepath, sectorString, str(camera))
            
                # only execute curl-download if required
                if not os.path.isfile(local_tar_filepath):
                    os.system(command)

                    # untar the file
                    currentDir = os.getcwd()
                    os.chdir(sector_directory)
                    command = 'tar xvf %s' %    local_tar_filepath
                    os.system(command)
                    os.chdir(currentDir)
    
            # get the names of all the downloaded fits files
            prf_filepaths.extend( glob.glob( sector_directory + '/*/*.fits') )

        # pass back a list of all of the downloaded prf files
        return prf_filepaths
    
    # -------------------------------------------------------------------------------------
    # Data directories / Storage / etc
    # -------------------------------------------------------------------------------------
    def _fetch_tess_data_directory(self):
        '''
            Create a sub-directory within self.local_dir
            Falls back to using self.local_dir
        '''
        
        # Define a subdirectory within the main data-directory
        tess_filepath = os.path.join(self._fetch_data_directory(), 'tess')
        
        # Attempt to create the desired sub-directory
        # Fall back to main data-directory if any error
        if not os.path.isdir(tess_filepath):
            try:
                os.mkdir(tess_filepath)
            except (Exception, OSError) as error :
                print(error)
                print('setting tess download directory to be same as parent directory: %r ' % self.local_dir )
                tess_filepath = self.local_dir
        return tess_filepath
    
    def _fetch_tess_prf_directory(self):
        '''
            Create a sub-directory within self._fetch_tess_data_directory()
            Falls back to using self._fetch_tess_data_directory()
            '''
        
        # Define a subdirectory within the main tess-directory
        prf_filepath = os.path.join(self._fetch_tess_data_directory(), 'prf')
        
        # Attempt to create the desired sub-directory
        # Fall back to main data-directory if any error
        if not os.path.isdir(prf_filepath):
            try:
                os.mkdir(prf_filepath)
            except (Exception, OSError) as error :
                print(error)
                print('setting tess-prf download directory to be same as parent directory: %r ' % self._fetch_tess_data_directory() )
                prf_filepath = self._fetch_data_directory()
        return prf_filepath

    
    def _ensure_tess_subdirectory_structure_exists(self, sectorNumber, cameraNumber, chipNumber):
        '''
            Ensure subdirectory structure exists
             - Subdirectories organized as tess_filepath/sector/camera/chip
        '''
        # Get the filepath to the tess data directory
        filepath = self.tess_dir
        
        # Attempt to recursively create sector/camera/chip subdirectories
        # Exit if these cannot be created
        for sub in [sectorNumber, cameraNumber, chipNumber]:
            filepath = os.path.join(filepath, sub )
            if not os.path.isdir(filepath):
                try:
                    os.mkdir(filepath)
                except (Exception, OSError) as error :
                    print(error)
                    sys.exit('Problem with sub-directory creation in ..._ensure_tess_subdirectory_structure_exists()...')
        
        return True

    

    def _fetch_tess_fits_filepath(self, fits_filename):
        '''
            Parses tess fits-filename and decides on destination filepath
            Filepath is just subdirectories within self.tess_dir
            Subdirectories organized as self.tess_dir/sector/camera/chip
            
            Input:
            ------
            fits_filename: name of fit file
             - looks like tess2018292175940-s0004-1-2-0124-s_ffic.fits
             
            Returns:
            --------
            dictionary 
             - looks like
             {
             'sectorNumber'  : 4,
             'cameraNumber'  : 1,
             'chipNumber'    : 2,
             'filepath'      : self.tess_dir/4/1/2/tess2018292175940-s0004-1-2-0124-s_ffic.fits
             }

             
        '''
        # can split on "-" to extract sector number, camera & chip
        sectorNumber   = int(fits_filename.split("-")[1][-4:])
        cameraNumber   = int(fits_filename.split("-")[2])
        chipNumber     = int(fits_filename.split("-")[3])
    
        # derive filepath
        subDir      = os.path.join(self.tess_dir, str(sectorNumber), str(cameraNumber), str(chipNumber) )
        filepath    = os.path.join(subDir, fits_filename)
        
        # ensure directory-structure exists
        self._ensure_tess_subdirectory_structure_exists(str(sectorNumber), str(cameraNumber), str(chipNumber))
        
        return {
            'sectorNumber'  : sectorNumber,
            'cameraNumber'  : cameraNumber,
            'chipNumber'    : chipNumber,
            'filepath'      : filepath
        }

    def _get_download_script(self , sectorNumber ):
        '''
            There are scripts available online to facilitate bulk downloads
            E.g. https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_4_ffic.sh
            
            This function will get the script that is relevant for a specific TESS sector
        '''
        filename    = 'tesscurl_sector_%s_ffic.sh' % sectorNumber
        filepath    = os.path.join(self.tess_dir , filename)
        if not os.path.isdir(filepath):
            try:
                os.system("curl -C - -L -o %s https://archive.stsci.edu/missions/tess/download_scripts/sector/%s" % ( filepath , filename))
            except (Exception, OSError) as error :
                print(error)
                sys.exit('Problem obtaining download script: %r' % filename )
        return filepath
    
    def _prf_filenames(self, cameraNumber, chipNumber):
        '''
            ...
        '''
        for sectorString in ['start_s0001','start_s0004']:
            for i in range(4):
                'https://archive.stsci.edu/missions/tess/models/prf_fitsfiles/%s/tess_prf_camera_%s.tgz' % (sectorString, str(camera))
    
    # -------------------------------------------------------------------------------------
    # curl commands for download
    # -------------------------------------------------------------------------------------

    @functools.lru_cache(maxsize=10)
    def _parse_download_script(self , shFile ):
        '''
            Parse a script downloaded using '_download_download_script'
             - I.e. read it to understand the names of the files that can be downloaded and hence where we want to file them away
            
            Contents look like ~16,000 rows of ...
            curl -C - -L -o tess2018292175940-s0004-1-2-0124-s_ffic.fits https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/tess2018292175940-s0004-1-2-0124-s_ffic.fits
            
            Will output a dictionary containing ...
             - list of fits files
            E.g. [ 'tess2018292095940-s0004-1-1-0124-s_ffic.fits' , '']
            
             - Will also parse the fits files to understand their sector, camera, chip, ...
            E.g. [(4,1,1,)]
            
             - Will also regenerate a 'corrected' curlcommand to put the file into the location we want 
            E.g. 'curl -C - -L -o ~/.shifty_data/tess/4/4/4/tess2018292105940-s0004-4-4-0124-s_ffic.fits https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/tess2018292105940-s0004-4-4-0124-s_ffic.fits'
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
    # The methods below are for the downloading of a limited, pre-defined set of TESS data
    # Intended for testing & development
    # -------------------------------------------------------------------------------------
    
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





class HSTDownloader(Downloader):
    '''
        You know we'll want to do it at some point !!!
    '''

    def __init__(self, ) :
            
        # - Allow ourselves to use any Downloader methods
        Downloader.__init__(self, )

        # - Local directory for saving data
        self.hst_dir = self._fetch_hst_data_directory()


    # -------------------------------------------------------------------------------------
    # Public methods for data download
    # -------------------------------------------------------------------------------------



    # -------------------------------------------------------------------------------------
    # Data directories / Storage / etc
    # -------------------------------------------------------------------------------------
    def _fetch_hst_data_directory(self):
        '''
            Create a sub-directory within self.local_dir
            Falls back to using self.local_dir
        '''

        # Define a subdirectory within the main data-directory
        hst_filepath = os.path.join(self._fetch_data_directory(), 'hst')

        # Attempt to create the desired sub-directory
        # Fall back to main data-directory if any error
        if not os.path.isdir(hst_filepath):
            try:
                os.mkdir(hst_filepath)
            except (Exception, OSError) as error :
                print(error)
                print('setting hst download directory to be same as parent directory: %r ' % self.local_dir )
                hst_filepath = self.local_dir
        return hst_filepath

