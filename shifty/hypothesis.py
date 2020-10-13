'''
    Definition of shifty.py's "OrbitHypothesis" class which we use to keep track of the different
    shifts we attempt. 
'''

# -------------------------------------------------------------------------------------
# Third party imports
# -------------------------------------------------------------------------------------
import os, sys
from collections import OrderedDict
import numpy as np
import functools


# -------------------------------------------------------------------------------------
# Any local imports
# -------------------------------------------------------------------------------------
# N/A


# -------------------------------------------------------------------------------------
# Various class definitions for generating/evaluating hypotheses in shifty
# -------------------------------------------------------------------------------------

class OrbitHypothesis():
    '''
        Container for pixel-shifts and alpha-dot, beta-dot, gamma, gamma-dot
    '''

    def __init__(self, shift_x, shift_y , alpha_dot, beta_dot, gamma, gamma_dot = 0 ):
        self.shift_x = shift_x
        self.shift_y = shift_y

        self.alpha_dot  = alpha_dot
        self.beta_dot   = beta_dot
        self.gamma      = gamma
        
        self.gamma_dot  = gamma_dot





class OrbitHypothesisGenerator(  ):
    '''
        Inputs:
        -------
        (d_min=100, d_max=300) or (alpha, beta, gamma), 
        
        ImageDataSet
        
        Methods:
        --------
        get_alpha_beta_gamma()
        get_hypotheses() => list of `OrbitHypothesis`

    '''
    
    def __init__(self, ImageDataSet , P=None, T=None ):
        
        # If P & T are *NOT* passed in as parameters, derive from IDS
        self.P = self._set_default_P_from_IDS(ImageDataSet) if P == None else P
        self.T = self._set_default_T_from_IDS(ImageDataSet) if T == None else T
    
    # -------------------------------------------------------------------------------------
    # Internal Methods
    # -------------------------------------------------------------------------------------
    
    def _set_default_P_from_IDS(self, ImageDataSet):
        ''' arc-seconds '''
        return 21.
    
    def _set_default_T_from_IDS(self, ImageDataSet):
        ''' days '''
        return 27.
    
    def _get_theta_max(self, gamma) :
        ''' Max Bound Velocity : arcsec/day '''
        return 360*3600*gamma**1.5 / 365.25 * np.sqrt(2)

    def _get_velocity_resolution(self,):
        ''' Following Bernstein and setting this to be ~P/T '''
        return self.P / self.T
    
    def _get_default_rate_of_motion_arrays(self, gamma):
        '''
            Set the default search range as a func of gamma
            alpha-dot_min = beta-dot_min == 0
            alpha-dot_max = beta-dot_max == theta_max(gamma)
            
            Then the range-lists span zero-to-max, in units of the velocity-resolution
        '''
        return np.arange(0 ,  self._get_theta_max(gamma) , self._get_velocity_resolution() )
    
    def _get_required_gamma_arrays():
        ''' Need to think this through more ... Bernstein '''
        gamma_max = 0. # = 1.0 / infinity
        gamma_max = 1./25.
        #N_gamma = 0.5 * (gamma/gamma_max) * (P/1)**(-1) * (T/1)
        #N_gamma = max(N_gamma, 1)
        return np.arange( [gamma_max] )

    # -------------------------------------------------------------------------------------
    # Public Methods
    # -------------------------------------------------------------------------------------
    
    def generate_hypotheses_over_gamma(self , gamma_list ):
        '''
            Will loop through supplied list of gamma-coefficients
            Will enumerate all possible alpha-dot, beta-dot values that will be investigated
            For each individual alpha-dot, beta-dot, gamma, it will generate the appropriate list of "shifts"
        '''

        hypotheses = []
        
        for gamma in gamma_list:
            
            # Get the min & max allowed rates of motion for a given gamma (distance)
            rate_of_motion_array = self._get_default_rate_of_motion_arrays(gamma)
        
            for alpha_dot in rate_of_motion_array:
                for beta_dot in rate_of_motion_array:
    
                    # Do something to generate shifts (in x & y pixels) as a function of alpha, beta, gamma
                    # - will presumably use ImageDataSet
                    
                    # Pass shifts to OrbitHypothesis object ...
                    hypotheses.append( OrbitHypothesis(shift_x, shift_y , alpha_dot, beta_dot, gamma) )


        return hypotheses
