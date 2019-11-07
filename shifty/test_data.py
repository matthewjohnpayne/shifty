"""
tests for shifty

to date, very few tests exist ...
"""


# -------------------------------------------------------------------------------------
# Third party imports
# -------------------------------------------------------------------------------------
import os, sys
import astropy
import numpy as np

# -------------------------------------------------------------------------------------
# Any local imports
# -------------------------------------------------------------------------------------
import data


# -------------------------------------------------------------------------------------
# Test "data" module
# -------------------------------------------------------------------------------------

def test_ImageLoader():
    ''' Test the ImageLoader parent class '''
    
    # Test creation of IL object
    IL = data.ImageLoader()
    assert isinstance(IL , data.ImageLoader), 'IL did not get created as expected'

    # Check has expected attributes
    assert 'local_dir' in IL.__dict__, ' local_dir not defined in IL'

    # Check save directory
    assert IL.local_dir == os.path.join(os.path.expanduser('~'), '.shifty_data'), 'local save dir not as expected'
    
    # Check method to load single image
    # - for now just checking that it returns an opened fits-file (formated as an hdulist)
    # - later we will want to perform more checks
    fits_filepath = os.path.join(IL.local_dir, str(4), str(1), str(1) , 'tess2018292095940-s0004-1-1-0124-s_ffic.fits')
    hdulist = IL._load_image(fits_filepath)
    assert isinstance(hdulist, astropy.io.fits.hdu.hdulist.HDUList), 'did not return expected object type'
    for key in ['PRIMARY']:
        assert key in hdulist, '%r not in hdulist' % key



    print('test_ImageLoader')



def test_TESSImageLoader():
    ''' Test the TESSImageLoader child class '''

    # Test creation of T object
    T = data.TESSImageLoader()
    assert isinstance(T , data.TESSImageLoader), 'T did not get created as expected'
    
    # Check T object has expected attributes
    assert 'local_dir' in T.__dict__, ' local_dir not defined in T'


    # ---- test methods related to defining storage directories -------------------------------------

    # Test method to derive tess_subdirectory_structure
    expectedDirectory = os.path.join(T.local_dir, str(4), str(4), str(4) )
    assert T._ensure_tess_subdirectory_structure_exists(str(4), str(4), str(4))
    assert os.path.isdir(expectedDirectory)
    
    # Test method to derive a directory filepath from a fits-filename
    fits_filename = 'tess2018292095940-s0004-1-1-0124-s_ffic.fits'
    result = T._fetch_tess_fits_filepath(fits_filename)
    assert isinstance(result, dict)
    for key in ['sectorNumber','cameraNumber','chipNumber', 'filepath' ]:
        assert key in result, \
            '%r not in result ' %  key
    expectedFilepath = os.path.join(T.local_dir, str(4), str(1), str(1) , fits_filename)
    assert result['filepath'] == expectedFilepath, \
        'filepath, %r, does not equal expectedFilepath, %r' % (result['filepath'], expectedFilepath)
    
    
    # Test methods to get download script(s)
    sectorNumber = 4
    expectedScriptAddress = os.path.join( T.local_dir , 'tesscurl_sector_%s_ffic.sh' % sectorNumber )
    '''
    if os.path.isfile(expectedScriptAddress): os.remove(expectedScriptAddress)
    assert not os.path.isfile(expectedScriptAddress), '%r did not get removed' % expectedScriptAddress
    outFilepath = T._download_download_script(sectorNumber)
    assert  os.path.isfile(expectedScriptAddress), '%r did not get downloaded' % expectedScriptAddress
    '''
    outFilepath = expectedScriptAddress

    # Test method to parse a downloaded script
    dataDict = T._parse_download_script(outFilepath)
    assert isinstance(dataDict, dict), 'not a dict '
    for key in ['fits_files','sectorNumbers','cameraNumbers','chipNumbers', 'filepaths' , 'curlCommands']:
        assert key in dataDict, '%r not in dataDict ' %  key

        if key == 'sectorNumbers':
            assert np.all( [ _ == sectorNumber for _ in dataDict[key] ] ), \
                'sectorNumbers are wrong'
        if key == 'cameraNumbers' or key == 'chipNumbers':
            assert np.all( [ _ in [1,2,3,4] for _ in dataDict[key] ] ), \
                'cameraNumbers/chipNumbers are wrong'
        if key in 'curlCommands':
            pass



    # ----------- test methods for getting some small amount of sample/test data -----------------------
    
    # test whether "_define_test_data()" returns a list of filepaths
    result = T._define_test_data()
    expectedLength = 10
    assert len(result) == expectedLength and np.all( [ T.local_dir in _ and ".fits" in _ for _ in result] ), \
        '_define_test_data() returns strange results: %r ' % result
    
    # test whether "_ensure_test_data_available_locally" actually ensures the required test data are available
    T._ensure_test_data_available_locally()
    assert np.all( [ os.path.isfile(_) for _ in T._define_test_data() ] ), \
        ' not all filepaths exist ...'
    
    # test the "_load_images()" method when supplied with 'DEV' option
    # - should return list pf HDUs
    result = T._load_images('DEV')
    assert isinstance(result, list) and len(result) == len(T._define_test_data()), \
        'returned result from T._load_test_data() not as expected'
    assert np.all( [isinstance(_, astropy.io.fits.hdu.hdulist.HDUList) for _ in result]), \
        'did not return expected object types in list'

    # test the *get_image_data_set()* method on the test-data set
    # - if we specify no cleaning functions, should be simple
    # - should generate an ImageDataSet object
    result = T.get_image_data_set( file_spec_container='DEV' )
    assert isinstance(result, data.ImageDataSet ), \
        '*get_image_data_set()* did not return an ImageDataSet'
    assert np.all( [ _ in result.__dict__ for _ in ['images', 'obs_code'] ] ),
        'ImageDataSet object did not have the expected attributes'



    print('\t test_TESSImageLoader passed ')




# -------------------------------------------------------------------------------------
# Test "hypothesis" module
# -------------------------------------------------------------------------------------








# Won't need these calls if use pytest/similar
#test_ImageLoader()
test_TESSImageLoader()

