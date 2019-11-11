'''
    Tests of Classes / Methods for downloading fits files (e.g. from MAST) to local disk
'''


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
import downloader


# -------------------------------------------------------------------------------------
# Test "downloader" module
# -------------------------------------------------------------------------------------



def test_Downloader():
    ''' Test the Downloader class '''
    print('\nWorking on test_Downloader() ...')

    # test basic creation of Downloader
    # -----------------------------------------

    # Test creation of Downloader object
    D = downloader.Downloader()
    assert isinstance(D , downloader.Downloader ), \
        'Downloader did not get created as expected'
    
    # Check D object has expected attributes
    assert 'local_dir' in D.__dict__, ' local_dir not defined in Downloader'


    # test methods related to defining storage directories
    # -----------------------------------------

    # we expect the local directory to be in '~/.shifty_data'
    expectedDirectory = os.path.join(os.path.expanduser('~'), '.shifty_data')
    assert D.local_dir == expectedDirectory, \
        'expectedDirectory [%r] != D.local_dir [%r]' % (expectedDirectory,D.local_dir)

    print('\t completed tests of test_Downloader ')



def test_TESSDownloader():
    ''' Test the TESSDownloader class '''
    print('\nWorking on test_TESSDownloader() ...')


    # test basic creation of Downloader
    # -----------------------------------------

    # Test creation of TESSDownloader object
    T = downloader.TESSDownloader()
    assert isinstance(T , downloader.TESSDownloader ), \
        'Downloader did not get created as expected'

    # Check T object has expected attributes
    assert 'local_dir' in T.__dict__, \
        ' local_dir not defined in TESSDownloader'
    assert 'tess_dir' in T.__dict__, \
        ' tess_dir not defined in TESSDownloader'



    # test methods related to defining storage directories
    # -----------------------------------------

    # we expect the local tess directory to be in '~/.shifty_data/tess'
    expectedDirectory = os.path.join(T._fetch_data_directory(), 'tess')
    assert T.tess_dir == T._fetch_tess_data_directory() == expectedDirectory, \
        'expectedDirectory [%r] != T.tess_dir [%r]' % (expectedDirectory, T.tess_dir)

    # Test method to derive tess_subdirectory_structure
    expectedDirectory = os.path.join(T.tess_dir, str(4), str(4), str(4) )
    assert T._ensure_tess_subdirectory_structure_exists(str(4), str(4), str(4))
    assert os.path.isdir(expectedDirectory)
    
    # Test method to derive a directory filepath from a fits-filename
    fits_filename   = 'tess2018292095940-s0004-1-1-0124-s_ffic.fits'
    result          = T._fetch_tess_fits_filepath(fits_filename)
    assert isinstance(result, dict)
    for key in ['sectorNumber','cameraNumber','chipNumber', 'filepath' ]:
        assert key in result, \
            '%r not in result ' %  key
    expectedFilepath = os.path.join(T.local_dir,
                                    'tess',
                                    str(4),
                                    str(1),
                                    str(1) ,
                                    fits_filename)
    assert result['filepath'] == expectedFilepath, \
        'filepath, %r, does not equal expectedFilepath, %r' % (result['filepath'], expectedFilepath)
    


    # test methods related to getting/parsing download scripts
    # -----------------------------------------

    # test methods to get download script(s)
    sectorNumber = 4
    expectedScriptAddress = os.path.join(T.tess_dir ,'tesscurl_sector_%s_ffic.sh' % sectorNumber )
    if False:
        if os.path.isfile(expectedScriptAddress): os.remove(expectedScriptAddress)
        assert not os.path.isfile(expectedScriptAddress), '%r did not get removed' % expectedScriptAddress
        outFilepath = T._get_download_script(sectorNumber)
        assert  os.path.isfile(expectedScriptAddress), '%r did not get downloaded' % expectedScriptAddress
    else:
        outFilepath = expectedScriptAddress
        
                                         
    # Test method to parse a downloaded script
    dataDict = T._parse_download_script(outFilepath)
    assert isinstance(dataDict, dict), 'not a dict '
    for key in ['fits_files','sectorNumbers','cameraNumbers','chipNumbers','filepaths' ,'curlCommands']:
        
        # check all keywords present
        assert key in dataDict, \
            '%r not in dataDict ' %  key
        # check sector-Number is as expected
        if key == 'sectorNumbers':
            assert np.all( [ _ == sectorNumber for _ in dataDict[key] ] ), \
                'sectorNumbers are wrong'
        # check camera/chip numbers are as expected
        if key == 'cameraNumbers' or key == 'chipNumbers':
            assert np.all( [ _ in [1,2,3,4] for _ in dataDict[key] ] ), \
                'cameraNumbers/chipNumbers are wrong'
        # check curl command is as expected
        if key in 'curlCommands':
            assert np.all( [ 'curl -C - -L -o ' in _ for _ in dataDict[key] ] ), \
                'curl commands are wrong'



    # test methods for getting sample/test data
    # -----------------------------------------
    
    # test whether "_define_test_data()" returns a list of filepaths
    result = T._define_test_data()
    expectedLength = 10
    assert len(result) == expectedLength and np.all( [ T.local_dir in _ and ".fits" in _ for _ in result] ), \
        '_define_test_data() returns strange results: %r ' % result

    # test *_download_test_data()* method

    # test whether "_ensure_test_data_available_locally" ...
    #... actually ensures the required test data are available
    T._ensure_test_data_available_locally()
    assert np.all( [ os.path.isfile(_) for _ in T._define_test_data() ] ), \
        ' not all filepaths exist ...'


    print('\t completed tests of test_TESSDownloader ')




    # test methods for getting big chunks of ffi data
    # -----------------------------------------
    sectorNumber, cameraNumber, chipNumber = 4,3,2
    expectedDirectory = os.path.join(T.tess_dir, str(sectorNumber), str(cameraNumber), str(chipNumber) )
    print('this download will take a while ...')
    T.download_chip(sectorNumber, cameraNumber, chipNumber)
    assert os.path.isdir(expectedDirectory), \
        'expectedDirectory [%r] does not exist' % expectedDirectory
    downloaded_files=glob.glob( os.path.join(expectedDirectory , '*.fits'))
    assert len(downloaded_files) > 1000, \
        'Insufficient files [%d] in %r ' % (len(downloaded_files), expectedDirectory )

def test_HSTDownloader():
    ''' Test the HSTDownloader class '''
    print('\nWorking on test_HSTDownloader() ...')
    
    
    # test basic creation of HSTDownloader
    # -----------------------------------------
    
    # Test creation of HSTDownloader object
    H = downloader.HSTDownloader()
    assert isinstance(H , downloader.HSTDownloader ), \
        'Downloader did not get created as expected'

    # Check H object has expected attributes
    assert 'local_dir' in H.__dict__, \
        ' local_dir not defined in HSTDownloader'
    assert 'hst_dir' in H.__dict__, \
        ' hst_dir not defined in HSTDownloader'



    # test methods related to defining storage directories
    # -----------------------------------------

    # we expect the local hst directory to be in '~/.shifty_data/hst'
    expectedDirectory = os.path.join(H._fetch_data_directory(), 'hst')
    assert H.hst_dir == H._fetch_hst_data_directory() == expectedDirectory, \
        'expectedDirectory [%r] != H.hst_dir [%r]' % (expectedDirectory, H.hst_dir)


    print('\t completed tests of test_HSTDownloader ')





# Won't need these calls if use pytest/similar
test_Downloader()
test_TESSDownloader()
test_HSTDownloader()

