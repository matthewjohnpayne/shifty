

# -------------------------------------------------------------------------------------
# Third party imports
# -------------------------------------------------------------------------------------
import os, sys
import astropy
from astropy.io import fits

# -------------------------------------------------------------------------------------
# Any local imports
# -------------------------------------------------------------------------------------
# N/A


# -------------------------------------------------------------------------------------
# Various class definitions for *Data Handling* in shifty
# -------------------------------------------------------------------------------------

class ImageDataSet():
    '''
        set of images w/ times & wcs which are suitable for stacking (e.g. stars have been masked, bad cadences removed, …)
        list of `ImageHDU` objects
        observatory_position (JPL code)
        get_observatory_barycentric_positions()
        get_theta_wcs()  expose the transformation of pixel to theta space

    '''

    def __init__(self, images , obs_code) :

        self.images = images
        self.obs_code = obs_code


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




class ImageLoader():
    '''
        Parent class for ...
         - TessImageLoader, HubbleImageLoader, PanstarrsImageLoader, ...
        
        => input: sector number, camera, chip, detector_range, bad_data_thresholds
        _load_images()
        _remove_stars()  (i.e. set pixels that contain stars to NaN)
        _remove_bad_cadences()
        _remove_scattered_light_problem_areas()
        _remove_strap_regions()
        get_image_data_set() => ImageDataSet
        

    '''
    
    def __init__(self, ) :
        pass

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

        data_dir = os.path.join(os.path.expanduser('~'), 'shifty_data')
        
        # If it doesn't exist, make a new data directory
        if not os.path.isdir(data_dir):
            
            try:
                os.makedirs(data_dir)
            
            # downloads locally if OS error occurs
            except OSError:
                warnings.warn('Warning: unable to create {}. '
                              'Download directory set to be the current working directory instead.'.format(data_dir))
                data_dir = '.'
                                                
        return data_dir

    # -------------------------------------------------------------------------------------
    # The methods below are for the loading of *general* fits-format data files
    # -------------------------------------------------------------------------------------

    def _load_images():
        '''
            Overall image loader as envisaged by Geert & Matt over coffee
        '''
        
    def _load_image(self , fits_file):
        '''
            Load a single image
            Basically a wrapper around astropy.fits.open
            But with the potential for added functionality (to be added later)
        '''
        return fits.open(os.path.join(self.local_dir, fits_file))



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
        
    def __init__(self, local_dir = None) :
        
        # - Allow ourselves to use ImageLoader methods
        super().__init__()
    
        # - Local directory for saving data
        self.local_dir = self._fetch_data_directory() if local_dir == None else local_dir

    # -------------------------------------------------------------------------------------
    # The methods below are for the "CLEANING" of TESS data
    # -------------------------------------------------------------------------------------
    def _mask_stars(self,):
        ''' 
            We want to remove stars in some way
            Barentsen & Payne discussed a simple mask: i.e. set pixels that contain stars to NaN
            This would be done based on GAIA positions
            
            This is *NOT* subtraction (see _subtract_stars below )
            
            Presumably only one of _mask_stars / _subtract_stars is required, but I am 100% certain that Holman will at least want to experiment with subtraction
            
        '''
        pass
    
    def _subtract_stars(self,):
        '''
            We want to remove stars in some way
            Holman & Payne have generally assumed some form of subtraction

            This is *NOT* masking (see _mask_stars )
            
            Presumably only one of _mask_stars / _subtract_stars is required, but I am 100% certain that Holman will at least want to experiment with subtraction

        '''
        pass

    def _remove_bad_cadences(self,):
        '''
            In many cases it may be most straightforward to simply eliminate entire exposures
            E.g. If they have terribly high stray-light across the entire exposure
        '''
        pass

    def _remove_scattered_light_problem_areas(self,):
        ''' 
            TESS has certain regions (of the chip(s)) in certain exposures that are known to have ...
            ... high levels of polluting scattered light
            We may want to completely mask these out
        '''
        pass

    def _remove_strap_regions(self,):
        '''
            TESS has certain regions (of the chip(s)) in which the backing-straps provide confusingly high signals
            We may want to completely mask these out (or perhaps flag them in some alternative manner)
        '''
        pass


    # -------------------------------------------------------------------------------------
    # The methods below are for the loading of a limited, pre-defined set of TESS data
    # Intended for the facilitation of testing
    # -------------------------------------------------------------------------------------

    def _load_test_data(self,):
        '''
            Load limited set of fits files for testing from local disk
        '''
        # If the test data doesn't exist on local disk, go and download the data from MAST
        self._ensure_test_data_available_locally()
    
        # Now that we are guaranteed that the data exists, we can open it
        # *** MJP: I Expect that the subsequent command/call will be updated (at some time in the future) to use a more general load function
        pass
    

    def _ensure_test_data_available_locally(self,):
        '''
            Check whether local versions of the fits files are available on local disk
            If they do not exist, download from MAST
        '''
        for fits_file in self._define_test_data():
            if not os.path.isfile( os.path.join(self.local_dir, fits_file) ):
                self._download_test_data(fits_file)
        print("\n\t All test data available locally ")

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
        return [
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







