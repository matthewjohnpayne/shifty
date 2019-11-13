'''
    methods to propagate known orbits for shifty.py
     - can be known *objects* or specified *orbit*
    
'''

# -------------------------------------------------------------------------------------
# Third party imports
# -------------------------------------------------------------------------------------
import os, sys
import wget
import numpy as np
import functools
import subprocess
import pandas as pd
import io
import csv


# -------------------------------------------------------------------------------------
# Any local imports
# -------------------------------------------------------------------------------------
from downloader import Downloader

# -------------------------------------------------------------------------------------
# Various class definitions for *propagating orbits* in shifty
# -------------------------------------------------------------------------------------

class Known(Downloader):
    '''
        class to deal with propagatation of known orbits
        
        methods:
        --------
        get_known_RADEC()
        convert_Keplerian_to_Alpha()
        convert_Alpha_to_Keplerian()
        _get_object_RADEC_from_horizons()
        _get_object_XYZ_from_horizons
        _get_orbit_RADEC()
        _convert_XYZ_to_RADEC()
        _propagate_orbit_keplerian()
        _propagate_orbit_nbody()
        
        
    '''
    
    def __init__(self, ) :
        pass
    
    # -------------------------------------------------------------------------------------
    # Public Methods
    # -------------------------------------------------------------------------------------
    def get_known_RADEC(self, **kwargs ):
        '''
           general method to get RADEC for specified object/orbit
           - calls specific sub-method as necessary
        '''
        if 'obs_code' not in kwargs or 'times' not in kwargs:
            sys.exit( 'get_known_RADEC() always requires at least *obs_code* and an array of *times* to be specified' )
        
        if 'object' in kwargs:
            self_get_object_RADEC_from_horizons(**kwargs)
        elif 'orbit' in kwargs:
            self._get_orbit_RADEC(**kwargs)
        else:
            sys.exit( 'do not know how to proceed to get RADEC' )
        
    def convert_Keplerian_to_Alpha(self, **kwargs ):
        '''
            convert a Keplerian orbit representation to the "tangent-plane" Alpha,Beta,Gamma,... representation
        '''
        pass
    def convert_Alpha_to_Keplerian(self, **kwargs ):
        '''
            inverse of convert_Keplerian_to_Alpha
        '''
        pass

    # -------------------------------------------------------------------------------------
    # Private Methods
    # -------------------------------------------------------------------------------------

    def _get_object_RADEC_from_horizons(self, objectName, obsCode, times ):
        '''
            query horizons for the RA,DEC of a *known object* at a sequence of times
            
        '''
        pass
    
    def _get_object_XYZ_from_horizons(self, objectName, obsCode, times):
        '''
            query horizons for the BARYCENTRIC POSITION of a *known object* at a sequence of times
        '''
    

    def _get_orbit_RADEC(self, **kwargs ):
        '''
            propagate a specified orbit and then calculate the expected RA,DEC at a sequence of times
             - wrapper around lower level method
            [[** to make life easier, demand observatory position be directly supplied ?? **]]
            
            inputs:
            -------
            
            returns:
            --------
        '''
        # do the propagation of the orbit
        if 'keplerian' in kwargs and kwargs['keplerian']:
            XYZs_ = self._propagate_orbit_keplerian()
        elif 'nbody' in kwargs and kwargs['nbody']:
            XYZs_ = self._propagate_orbit_nbody()
        else:
            sys.exit('do not know how to propagate ... ')

        # convert the orbit to apparent RA, Dec
        return self._convert_XYZ_to_RADEC(XYZs_)


    def _propagate_orbit_keplerian():
        '''
            propagate an orbit using keplerian methods
            [[think about chebyshev]]
        '''
        pass
    def _propagate_orbit_nbody():
        '''
            propagate an orbit using nbody methods
            [[think about chebyshev]]
        '''
        pass
    def _convert_XYZ_to_RADEC():
        ''' 
            convert XYZ positions to apparent RA, Dec
        '''
        pass
