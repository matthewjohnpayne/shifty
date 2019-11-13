'''
    methods to deal with pixel-related operations in shify.py
    includes 
     - WCS-related transformations : pixel <-> (ra,dec )
     - pixel-adjacency lookups/searches ...
'''

# -------------------------------------------------------------------------------------
# Third party imports
# -------------------------------------------------------------------------------------
import os, sys
import numpy as np


import astropy
from astropy.time import Time
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord


# -------------------------------------------------------------------------------------
# Any local imports
# -------------------------------------------------------------------------------------
from downloader import Downloader

# -------------------------------------------------------------------------------------
# Various class definitions for *pixel handling* in shifty
# -------------------------------------------------------------------------------------

class Pixel(Downloader):
    '''
        class to deal with to deal with pixel-related operations
        
        methods:
        --------
        get_per_pixel_RADEC()
        adjacent_pixels()
        radius_search_over_pixels()
        
    '''
    
    def __init__(self, ) :
        
        # - Allow ourselves to use any XXX methods
        # XXX.__init__(self, )
        pass

    # -------------------------------------------------------------------------------------
    # Public Methods
    # -------------------------------------------------------------------------------------
    def get_per_pixel_RADEC(self,  header, data ):
        '''
            use wcs to get RADEC of each pixel
            https://docs.astropy.org/en/stable/api/astropy.wcs.WCS.html#astropy.wcs.WCS.all_pix2world
            
            NOTE: axis ordering difference between numpy array and FITS
            
            inputs:
            -------
            header
            - header portion of fits, must contain wcs solution
            data
            - data portion of fits: np.ndarray
            
            returns:
            --------
            astropy.coordinates.sky_coordinate.SkyCoord object
            [[could add option to return 2 separate ra,dec arrays ]]
        '''
        # not sure whether ultimately necessary, but for now enforce 2D-array
        assert len(data.shape)==2, 'data is not 2D:: data.shape=%r' % data.shape
        
        # set up a meshgrid of pixel addresses for the image-data pixels
        xx, yy = np.meshgrid(range(np.shape(data)[1]), range(np.shape(data)[0]) )
        
        # calculate the RA,DEC values for each pixel using astropy-wcs
        ra,dec = np.array(WCS(header).all_pix2world(xx, yy, 1))
        
        # return an astropy SkyCoord object
        return SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')

    def adjacent_pixels(self, ):
        '''
           given a pixel address, return the pixel addresses of the adjacent pixels 
           account for edge-effects  
           
           e.g.:
           -----
           
                x
           y    0       1       2       3
           -    --------------------------
           2   (0,2)    (1,2)   (2,2)   (3,2)
           1   (0,1)    (1,1)   (2,1)   (3,1)
           0   (0,0)    (1,0)   (2,0)   (3,0)
           
           (a): (1,1) has adjacencies [(0,0),(1,0),(2,0),(0,1),(2,1),(0,2),(1,2),(2,2)]
           (b): (0,0) has adjacencies [(1,0),(0,1),(1,1)]
           
           inputs:
           -------
           
           returns:
           --------
        '''

    def radius_search_over_pixels(self, ):
        '''
           given a precise ra,dec, return the pixels touched by a circle of radius, r
           
           inputs:
           -------
           
           returns:
           --------
        '''
        # map ra,dec to single pixel
        # get adjacent pixels
        # get ra, dec (or uv) of adjacent pixels
        # get mid-points between adjacent pixels (account for ra ~ 360)
        # get angular separation between precise ra,dec (from input) and mid-points
        # return pixels within radius, r
        pass

    # -------------------------------------------------------------------------------------
    # Private Methods
    # -------------------------------------------------------------------------------------
