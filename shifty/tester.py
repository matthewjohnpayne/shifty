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
    
    # Check has expected attributes
    assert 'local_dir' in T.__dict__, ' local_dir not defined in T'

    # Test method to derive tess_subdirectory_structure
    expectedDirectory = os.path.join(T.local_dir, str(4), str(4), str(4) )
    assert T._ensure_tess_subdirectory_structure_exists(str(4), str(4), str(4))
    assert os.path.isdir(expectedDirectory)
    
    # Test method to derive a directory filepath from a fits-filename
    fits_filename = 'tess2018292095940-s0004-1-1-0124-s_ffic.fits'
    result = T._fetch_tess_fits_filepath(fits_filename)
    assert isinstance(result, dict)
    for key in ['sectorNumber','cameraNumber','chipNumber', 'filepath' ]:
        assert key in result, '%r not in result ' %  key
    expectedFilepath = os.path.join(T.local_dir, str(4), str(1), str(1) , fits_filename)
    assert result['filepath'] == expectedFilepath, 'filepath, %r, does not equal expectedFilepath, %r' % (result['filepath'], expectedFilepath)
    
    
    # Test methods to get download script(s)
    sectorNumber = 4
    expectedScriptAddress = os.path.join( T.local_dir , 'tesscurl_sector_%s_ffic.sh' % sectorNumber )
    print(expectedScriptAddress)
    if os.path.isfile(expectedScriptAddress): os.remove(expectedScriptAddress)
    assert not os.path.isfile(expectedScriptAddress), '%r did not get removed' % expectedScriptAddress
    outFilepath = T._download_download_script(sectorNumber)
    assert  os.path.isfile(expectedScriptAddress), '%r did not get downloaded' % expectedScriptAddress

    # Test method to parse the download file
    dataDict = T._parse_download_script(outFilepath)
    assert isinstance(dataDict, dict), 'not a dict '
    for key in ['fits_files','sectorNumbers','cameraNumbers','chipNumbers', 'filepaths' , 'curlCommands']:
        assert key in dataDict, '%r not in dataDict ' %  key

        if key == 'sectorNumbers':
            assert np.all( [ _ == sectorNumber for _ in dataDict[key] ] ), 'sectorNumbers are wrong'
        if key == 'cameraNumbers' or key == 'chipNumbers':
            assert np.all( [ _ in [1,2,3,4] for _ in dataDict[key] ] ), 'cameraNumbers/chipNumbers are wrong'
        if key in 'curlCommands':
            print(dataDict[key][:2])

    print('test_TESSImageLoader')





# -------------------------------------------------------------------------------------
# Test "hypothesis" module
# -------------------------------------------------------------------------------------








# Won't need these calls if use pytest/similar
test_ImageLoader()
test_TESSImageLoader()

