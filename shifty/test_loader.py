"""
    Test classes / methods to load fits file from disk into an "ImageDataSet" object
"""


# -------------------------------------------------------------------------------------
# Third party imports
# -------------------------------------------------------------------------------------
import os, sys
import astropy
import numpy as np
from astropy.time import Time

# -------------------------------------------------------------------------------------
# Any local imports
# -------------------------------------------------------------------------------------
import loader
import data

# -------------------------------------------------------------------------------------
# Test "loader" module
# -------------------------------------------------------------------------------------



def test_ImageLoader():
    ''' Test the ImageLoader parent class '''
    print('\nWorking on test_ImageLoader() ...')


    # test basic creation of ImageLoader
    # -----------------------------------------

    # Test creation of IL object
    IL = loader.ImageLoader()
    assert isinstance(IL , loader.ImageLoader), 'IL did not get created as expected'

    # Check has expected attributes
    assert 'local_dir' in IL.__dict__, ' local_dir not defined in IL'

    # Check save directory
    assert IL.local_dir == os.path.join(os.path.expanduser('~'), '.shifty_data'), 'local save dir not as expected'



    # test image-loading method(s)
    # -----------------------------------------

    # Check method to load single image
    # - for now just checking that it returns an opened fits-file
    #   (formated as an hdulist)
    # - later we will want to perform more checks
    fits_filepath = os.path.join(IL.local_dir,
                                 'tess',
                                 str(4), # sector
                                 str(1), # camera
                                 str(1), # chip
                                 'tess2018292095940-s0004-1-1-0124-s_ffic.fits')
    hdulist = IL._load_image(fits_filepath)
    assert isinstance(hdulist, astropy.io.fits.hdu.hdulist.HDUList), \
        'did not return expected object type'
    for key in ['PRIMARY']:
        assert key in hdulist, '%r not in hdulist' % key



    print('\t completed tests of test_ImageLoader')





def test_TESSImageLoader():
    ''' Test the TESSImageLoader child class '''
    print('\nWorking on test_TESSImageLoader() ...')

    # test basic creation of TESSImageLoader
    # -----------------------------------------

    # Test creation of TESSImageLoader object
    T = loader.TESSImageLoader()
    assert isinstance(T , loader.TESSImageLoader), 'T did not get created as expected'
    
    # Check T object has expected attributes
    assert 'obs_code' in T.__dict__, ' obs_code not defined in T'






    # test methods for loading TESS data
    # -----------------------------------------
    
    # test the "_load_images()" method when supplied with 'DEV' option
    # - should return list pf HDUs
    result = T._load_images( development = True )
    assert isinstance(result, list) and len(result) == len(T._define_test_data()), \
        'returned result from T._load_test_data() not as expected'
    assert np.all( [isinstance(_, astropy.io.fits.hdu.hdulist.HDUList) for _ in result]), \
        'did not return expected object types in list'


    # test the convenience function *_load_test_images()*
    # this should be identical to using _load_images(  development = True )
    result = T._load_test_images()
    assert isinstance(result, list) and len(result) == len(T._define_test_data()), \
        'returned result from T._load_test_data() not as expected'
    assert np.all( [isinstance(_, astropy.io.fits.hdu.hdulist.HDUList) for _ in result]), \
        'did not return expected object types in list'

    # test the *_load_images()* method specifying sectorNumber/cameraNumber/chipNumber
    sectorNumber, cameraNumber, chipNumber = 4,3,2
    # result = T._load_images( sectorNumber=sectorNumber, cameraNumber=cameraNumber, chipNumber=chipNumber )
    print(' *NO* test of _load_images specifying sectorNumber/cameraNumber/chipNumber in place ...')

    # test the *_load_sector_camera_chip()* method
    # this should be identical to using _load_images()
    sectorNumber, cameraNumber, chipNumber = 4,3,2
    # result = T._load_sector_camera_chip( sectorNumber, cameraNumber, chipNumber)
    print(' *NO* test of _load_sector_camera_chip in place ...')




    # test methods for turning fits-file data into ImageDataSet objects
    # -----------------------------------------

    # test the *get_image_data_set()* method on the test-data set
    # - if we specify no cleaning functions, should be simple
    # - should generate an ImageDataSet object
    result = T.get_image_data_set( **{'development' : True} )
    assert isinstance(result, data.ImageDataSet ), \
        '*get_image_data_set()* did not return an ImageDataSet'
    assert np.all( [ _ in result.__dict__ for _ in ['headers', 'data','unc', 'obs_code'] ] ),\
        'ImageDataSet object did not have the expected attributes'


    # repeat the test of *get_image_data_set()* method on the test-data set
    # - but now pass a cleaning_dict with params all set False
    # - should have the same outcome as above
    cpd = { _ : False for _ in ['mask' ,'subtract', 'bad_cad', 'scat', 'strap' ] }
    result = T.get_image_data_set( **{'development' : True} , **cpd )
    assert isinstance(result, data.ImageDataSet ), \
        '*get_image_data_set()* did not return an ImageDataSet'
    assert np.all( [ _ in result.__dict__ for _ in ['headers', 'data','unc', 'obs_code'] ] ),\
        'ImageDataSet object did not have the expected attributes'
    



    # test methods for loading pixel-response-function data
    # -----------------------------------------
    print('prfs ...', T.download_prf() )
    result = T.get_prfs()
    assert isinstance(result, list), \
        'not a list: %r' % type(result)
    for item in result:
        assert isinstance(item, astropy.io.fits.hdu.hdulist.HDUList), \
            'not an HDU object: %r' % type(item)





    # test methods for masking data
    # -----------------------------------------
    #result = T.get_image_data_set( **{'development' : True , 'mask':True }  )



    print('\t completed tests of test_TESSImageLoader ')




# Won't need these calls if use pytest/similar
test_ImageLoader()
test_TESSImageLoader()

