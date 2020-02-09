"""
    Test classes / methods to to prepare fits files
"""


# -------------------------------------------------------------------------------------
# Third party imports
# -------------------------------------------------------------------------------------
import os, sys
import astropy
import numpy as np
import glob

import astropy
from astropy.io import fits
from astropy.time import Time

# -------------------------------------------------------------------------------------
# Any local imports
# -------------------------------------------------------------------------------------
sys.path.append(os.path.join(os.path.split(os.getcwd())[0], 'shifty'))
import preparer
import data

# -------------------------------------------------------------------------------------
# Test "preparer" module
# -------------------------------------------------------------------------------------



def test_ImagePreparer():
    ''' Test the ImagePreparer parent class '''
    print('\nWorking on test_ImagePreparer() ...')


    # test basic creation of ImagePreparer
    # -----------------------------------------

    # Test creation of IL object
    IL = preparer.ImagePreparer()
    assert isinstance(IL , preparer.ImagePreparer), 'IL did not get created as expected'

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



    print('\t completed tests of test_ImagePreparer')





def test_TESSImagePreparer():
    ''' Test the TESSImagePreparer child class '''
    print('\nWorking on test_TESSImagePreparer() ...')

    # test basic creation of TESSImagePreparer
    # -----------------------------------------

    # Test creation of TESSImagePreparer object
    T = preparer.TESSImagePreparer()
    assert isinstance(T , preparer.TESSImagePreparer), 'T did not get created as expected'

    # Check has expected inherited attributes
    assert 'local_dir' in T.__dict__, \
        ' local_dir not defined in IL'
    assert T.local_dir == os.path.join(os.path.expanduser('~'), '.shifty_data'), \
        'local save dir not as expected'

    # Check T object has its own expected attributes
    assert 'obs_code' in T.__dict__, \
        ' obs_code not defined in T'





    # test methods for loading TESS FFI data
    # -----------------------------------------

    # (i) test the private-method *_parse_filespec()*
    # (i)(a) should get simple pass-through of file-list when 'fits_filepaths' used
    result = T._parse_filespec( fits_filepaths = T._ensure_test_data_available_locally() )
    assert result == T._ensure_test_data_available_locally(), \
        'did not get simple pass-through of fits_filepaths'

    # (i)(b) should get filepaths for development files when 'development' = True is specified
    result = T._parse_filespec( development = True )
    assert result == T._ensure_test_data_available_locally(), \
        'did not get development fits_filepaths'

    # (i)(c) should get all filepaths in chip sub-dir when 'sectorNumber', 'cameraNumber', 'chipNumber' specified
    sectorNumber, cameraNumber, chipNumber = 4,3,2
    result = T._parse_filespec(sectorNumber = sectorNumber,
                               cameraNumber = cameraNumber,
                               chipNumber   = chipNumber )
    result.sort()
    expected_directory_path = os.path.join( T.tess_dir,str(sectorNumber),str(cameraNumber),str(chipNumber))
    expected_fits_filepaths = glob.glob( os.path.join(expected_directory_path , '*.fits') )
    expected_fits_filepaths.sort()
    assert result == expected_fits_filepaths, \
        'did not get expected_fits_filepaths for %r' % expected_directory_path



    # (ii) test the private-method *_parse_patchspec()*
    # (ii)(a) if use pythonic-zero-based-numbering, then numbers unaltered
    X0,X1,Y0,Y1 = 11,17,1,4
    x0,x1,y0,y1 = T._parse_patchspec(patch=True, python=True, xlim=(X0,X1), ylim=np.array([Y0,Y1]) )
    assert [x0,x1,y0,y1]==[X0,X1,Y0,Y1], \
        'returned limits are not as expected ... %r' % [x0,x1,y0,y1]

    # (ii)(b) if use pixel-one-based-numbering, then numbers shifted by one
    X0,X1,Y0,Y1 = 11,17,1,4
    x0,x1,y0,y1 = T._parse_patchspec(patch=True, pixel=True, xlim=(X0,X1), ylim=np.array([Y0,Y1]) )
    assert [x0,x1,y0,y1]==[X0-1,X1-1,Y0-1,Y1-1], \
        'returned limits are not as expected ... %r' % [x0,x1,y0,y1]



    # (iii) test the private-method *_load_test_images()*
    # - convenience function to load a small, pre-defined sample of test data
    # - should return list of HDUlists
    result = T._load_test_images()
    assert isinstance(result, list) and len(result) == len(T._define_test_data()), \
        'returned result from T._load_test_images() not as expected'
    assert np.all( [isinstance(_, astropy.io.fits.hdu.hdulist.HDUList) for _ in result]), \
        'did not return expected object types in list'


    # (iv) test the public-method *generate_cleaned_stack_file()*
    # - reads, cleans, saves to stack-file

    # (iv)(a) simplest use: make small stack file out of development-data, doing NO cleaning
    filepath = T.generate_cleaned_stack_file( development = True )
    assert isinstance(filepath, str) and os.path.isfile(filepath), \
        'returned value, %r, does not appear to be a valid filepath' % filepath
    assert isinstance( fits.open(filepath), astropy.io.fits.hdu.hdulist.HDUList), \
        'opening %r with astropy.fits.open does not work as expected' % filepath


    # (iv)(b) ... ***NEED MANY MORE TESTS HERE ***



    # (v) test the private-method *_generate_test_stack_file()*
    # - convenience function to load a small, pre-defined sample of test data into a single stack-file
    # - depends on generate_cleaned_stack_file() &  _parse_filespec()
    filepath = T._generate_test_stack_file()
    assert isinstance(filepath, str) and os.path.isfile(filepath), \
        'returned value, %r, does not appear to be a valid filepath' % filepath
    assert isinstance( fits.open(filepath), astropy.io.fits.hdu.hdulist.HDUList), \
        'opening %r with astropy.fits.open does not work as expected' % filepath




    # test methods for handling an overall "stack" fits-file
    # -----------------------------------------

    # (i) test the private-method *_initialize_stack_HDUlist()*
    # - should create a fits-file and an hdulist
    stack_filepath = os.path.join(T._fetch_data_directory(), 'stack.fits')
    if os.path.isfile(stack_filepath):
        os.remove(stack_filepath)
    assert not os.path.isfile(stack_filepath)
    result = T._initialize_stack_HDUlist(stack_filepath)
    assert isinstance( result, astropy.io.fits.hdu.hdulist.HDUList)
    assert os.path.isfile(stack_filepath)




    # test methods for handling TESS pixel-response-function data
    # -----------------------------------------

    # () test the public-method *get_prfs()*
    print(' *NO* test of get_prfs() in place ...')




    # test methods for parsing TESS FFI data
    # -----------------------------------------

    # (i) test the private-method *_get_midtime()*
    # - expect it to return an astropy.time.Time object with a single float value
    list_of_HDUs = T._load_test_images()
    header = list_of_HDUs[0][1].header
    result = T._get_midtime( header )
    assert isinstance( result , astropy.time.Time ), \
        'not an astropy Time object'
    assert isinstance( result.value , float ), \
        'not a float'
    expectedTime =      header['BJDREFI']    \
                +       header['BJDREFF']    \
                    +   0.5*(header['TSTART']      + header['TSTOP'])
    assert  result.value == expectedTime, \
        'not expected time value'



    # test methods for "cleaning" TESS FFI data
    # -----------------------------------------

    # (i) test the private-method *_clean_data()*
    # - Wrapper around the various "cleaning" methods below

    # (ii) test the private-method *_mask_stars()*

    # (ii)(a) demonstrate that the kwargs get updated in the desired manner
    kwargs = { 'dummy_key' : True , 'refcat_dict' : {} }
    list_of_HDUs = T._load_test_images()
    header  = list_of_HDUs[0][1].header
    data    = list_of_HDUs[0][1].data
    print(0, kwargs)
    T._mask_stars(header , data, **kwargs)
    print(1, kwargs)



    # () test the private-method *_subtract_stars()*

    # () test the private-method *_remove_bad_cadences()*

    # () test the private-method *_remove_scattered_light_problem_areas()*

    # () test the private-method *_remove_strap_regions()*





    print('\t completed tests of test_TESSImagePreparer ')
    print('\t\t N.B. tests are *INCOMPLETE*')




# Won't need these calls if use pytest/similar
#test_ImagePreparer()
test_TESSImagePreparer()

