"""
    tests for refcat download within shifty
    
"""


# -------------------------------------------------------------------------------------
# Third party imports
# -------------------------------------------------------------------------------------
import os, sys
import numpy as np


# -------------------------------------------------------------------------------------
# Any local imports
# -------------------------------------------------------------------------------------
from refcat import RefCat
import loader
# -------------------------------------------------------------------------------------
# Test "refcat" classes/modules
# -------------------------------------------------------------------------------------



def test_fetch_refcat_data_directory():
    
    # Test creation of RefCat object
    # -----------------------------------------
    RC = RefCat()
    assert isinstance(RC , RefCat ), \
        'RefCat did not get created as expected'
    for k,v in RC.__dict__.items():print(k,v)
    # Check RefCat object has expected attributes
    # -----------------------------------------
    assert 'local_dir' in RC.__dict__, \
        ' local_dir not defined in RefCat'
    assert 'refcat_dir' in RC.__dict__, \
        ' refcat_dir not defined in RefCat'

    # We expect the local refcat directory to be in '~/.shifty_data/refcat'
    expectedDirectory = os.path.join(RC._fetch_data_directory(), 'refcat')
    assert RC.refcat_dir == RC._fetch_refcat_data_directory() == expectedDirectory, \
        'expectedDirectory [%r] != RC.tess_dir [%r]' % (expectedDirectory, RC.tess_dir)

    print('\t Successfully tested *_fetch_refcat_data_directory*')



def test_download_refcat():
    
    # Test creation of RefCat object
    # -----------------------------------------
    RC = RefCat()
    assert isinstance(RC , RefCat ), \
        'RefCat did not get created as expected'

    # run the download script
    # -----------------------------------------
    RC.download_refcat(16)

    # check that the expected file is present in the refcat directory
    expected_file = os.path.join(RC.refcat_dir ,
                                 'hlsp_atlas-refcat2_atlas_ccd_00-m-16_multi_v1_cat.tbz')
    assert os.path.isfile( expected_file ), \
        'expected file does not exist ...'


    print('\t Successfully tested *download_refcat()*')


def test_untar_refcat():
    
    # Test creation of RefCat object
    # -----------------------------------------
    RC = RefCat()
    assert isinstance(RC , RefCat ), \
        'RefCat did not get created as expected'
    
    # run the download script
    # -----------------------------------------
    RC._untar_refcat(16)

    # check that the expected file is present in the refcat directory
    expected_file = os.path.join(RC.refcat_dir ,
                             '00_m_16')
    assert os.path.isdir( expected_file ), \
        'expected directory does not exist ...'
    
    
    print('\t Successfully tested *_untar_refcat()*')




def test_download_refcode():
    
    # Test creation of RefCat object
    # -----------------------------------------
    RC = RefCat()
    assert isinstance(RC , RefCat ), \
        'RefCat did not get created as expected'
    
    # run the download script
    # -----------------------------------------
    RC.download_refcode()

    # check that the expected file is present in the refcat directory
    expected_file = os.path.join(RC.refcat_dir ,
                                 'refcat.c')
    assert os.path.isfile( expected_file ), \
        'expected file does not exist ...'
    
    
    print('\t Successfully tested *download_refcode()*')


def test_compile_refcat():
    
    # Test creation of RefCat object
    # -----------------------------------------
    RC = RefCat()
    assert isinstance(RC , RefCat ), \
        'RefCat did not get created as expected'
    
    # run the compile script
    # -----------------------------------------
    RC.compile_refcat()

    # check that the expected file is present in the refcat directory
    expected_file = os.path.join(RC.refcat_dir ,
                             'refcat')
    assert expected_file == RC.refcat_filepath, \
        'returned filepath [%r] differs from expected [%r]' % ( RC.refcat_filepath , expected_file )
    assert os.path.isfile( expected_file ), \
        'expected file does not exist ...'

    print('\t Successfully tested *compile_refcat()*')


def test_read_refcat():
    
    # Test creation of RefCat object
    # -----------------------------------------
    RC = RefCat()
    assert isinstance(RC , RefCat ), \
        'RefCat did not get created as expected'

    # run the _read_refcat script
    # -----------------------------------------
    # The following is taken from refcat.man
    # - Report all stars within a radius of 1 deg from RA, Dec 180,10.
    # - Request input file extension ".rc2" explicitly and a header line
    # - (1064 stars, output in ATLAS format):
    # >>> refcat 180 10 -rad 1 -dir 00_m_16 -exten rc2 -hdr
    #
    # Query returns results which look like ...
    # 179.755775  10.937336 11.690 11.243 11.109 11.059 10.142 11.472 11.182
    # -----------------------------------------
    ra = 180
    dec = 10
    rad = 1
    mlim = 17
    code = RC.refcat_filepath
    dir = os.path.join(RC.refcat_dir, '00_m_16')
    result = RC._read_refcat(ra, dec, code, dir,
                         rad=rad , mlim=mlim )
    assert isinstance(result, dict), '%r not a dict ' % type(result)
    assert len(result) == 1064

    print('\t Successfully tested *_read_refcat()*')


def test_find_all_stars_on_image():
    
    # Test creation of RefCat object
    # -----------------------------------------
    RC = RefCat()
    assert isinstance(RC , RefCat ), \
        'RefCat did not get created as expected'
    
    # Need to create header, image_data to pass into *find_all_stars_on_image()*
    # -----------------------------------------
    T = loader.TESSImageLoader()
    cpd = { _ : False for _ in ['mask' ,'subtract', 'bad_cad', 'scat', 'strap' ] }
    IDS = T.get_image_data_set( **{'development' : True} , **cpd )
    
    # run the find_all_stars_on_image script
    # -----------------------------------------
    RC.find_all_stars_on_image(IDS.headers[0] , IDS.data[0])

    print('\t Successfully tested *find_all_stars_on_image()*')


# Won;t need these if I bother to get pytest working ...
"""
test_fetch_refcat_data_directory()
test_download_refcat()
test_untar_refcat()
test_download_refcode()
test_compile_refcat()
"""
test_find_all_stars_on_image()

