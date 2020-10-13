'''
   Classes / methods to load fits files into an "ImageDataSet" object
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
import imagedata

# -------------------------------------------------------------------------------------
# Various class definitions for *data import * in shifty
# -------------------------------------------------------------------------------------

class Loader(Downloader):
    '''
        (1) Loads fits-files
        (2) ...
        
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






class TESSLoader(ImagePreparer, TESSDownloader):
    '''
        For the specific preparation of *TESS* data 
        
        (1) Loads fits-files
        (2)...
        
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
    # "Patch" Related
    # -------------------------------------------------------------------------------------
    

    def _patch_to_limits(self, patch=1, patchsize=512):
        """Returns the corners (xmin, xmax, ymin, ymax) of a patch."""
        patch_x = int(patch % np.floor(2048 / patchsize))
        patch_y = int(np.floor((patch * patchsize) / 2048))
        return (self.XRANGE[0] + patch_x * patchsize,
                self.XRANGE[0] + (patch_x + 1) * patchsize,
                self.YRANGE[0] + patch_y * patchsize,
                self.YRANGE[0] + (patch_y + 1) * patchsize)
    
    def get_patch(self, patch=1, patchsize=512):
        """Returns an `ImageDataSet` object containing the pixel data in a patch."""
        xmin, xmax, ymin, ymax = self._patch_to_limits(patch=patch, patchsize=patchsize)
        images = np.empty((len(self.filenames), ymax-ymin, xmax-xmin))
        headers = []
        for idx, fn in tqdm(enumerate(self.filenames), total=len(self.filenames), desc='Reading FFIs'):
            with fitsio.FITS(fn) as fts:
                images[idx] = fts[1][ymin:ymax, xmin:xmax]
                headers.append(fts[1].read_header())
    
        # Return an ImageDataSet object
        return imagedata.ImageDataSet(images, headers)

    def get_high_quality_patch(self, quality_level=1, patch=1, patchsize=512):
        """ 
            As per *get_patch()*, returns an `ImageDataSet` object containing the pixel data in a patch.
            However, get_high_quality_patch will also remove poor "quality" exposures according to various metrics
        """
        # 0 => Keep everything
        if quality_level == 0:
            pass
        # 1 => Keep only exposures with tasoc_quality ...
        elif quality_level == 1:
            pass
        # 2 => Keep only exposures with tasoc_quality & low per-pixel variability
        elif quality_level == 2:
            pass
        else
            print('Could not understand quality level: %r' % quality_level )
            

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
    # PRF/PSF Related
    # -------------------------------------------------------------------------------------


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
            
            


    # -------------------------------------------------------------------------------------
    # The methods below are for the "PARSING" of TESS HDUS from FITS
    # -------------------------------------------------------------------------------------

    # get midtimes
    #midtimes = np.array( [self._get_midtime(item.header)  for item in hdul  \
    #                  if isinstance(item, astropy.io.fits.hdu.image.ImageHDU )] )

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
    




class HSTLoader(ImagePreparer):
    '''
        You know we'll want to do it at some point !!!
        '''
    pass








