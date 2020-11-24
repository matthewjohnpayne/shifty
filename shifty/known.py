'''
    methods to propagate known orbits for shifty.py
     - can be known *objects* or specified *orbit*

'''

# -----------------------------------------------------------------------------
# Third party imports
# -----------------------------------------------------------------------------
import os
import sys
import numpy as np
from astroquery.jplhorizons import Horizons


# -----------------------------------------------------------------------------
# Any local imports
# -----------------------------------------------------------------------------
sys.path.append(os.path.dirname(os.path.dirname(
                os.path.realpath(__file__))))
from shifty.downloader import Downloader


# -----------------------------------------------------------------------------
# Constants and Test data
# -----------------------------------------------------------------------------

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(
                        os.path.realpath(__file__))), 'dev_data')


# -----------------------------------------------------------------------------
# Various class definitions for *propagating orbits* in shifty
# -----------------------------------------------------------------------------

class Known(Downloader):  # MA: Why is Downloader a base class here???
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

        _interpolate_radec_for_sedna()


    '''

    def __init__(self, ):
        super().__init__()

    # -------------------------------------------------------------------------
    # Public Methods
    # -------------------------------------------------------------------------
    def get_known_RADEC(self, **kwargs):
        '''
           General method to get RADEC for specified object/orbit
           - calls specific sub-method as necessary
        '''
        if 'obs_code' not in kwargs or 'times' not in kwargs:
            sys.exit('get_known_RADEC() always requires at least '
                     '*obs_code* and an array of *times* to be specified')

        if 'object' in kwargs:
            self._get_object_RADEC_from_horizons(**kwargs)
        elif 'orbit' in kwargs:
            self._get_orbit_RADEC(**kwargs)
        else:
            sys.exit('Do not know how to proceed to get RADEC')

    def convert_Keplerian_to_Alpha(self, **kwargs):
        '''
            Convert a Keplerian orbit representation to the "tangent-plane"
            Alpha,Beta,Gamma,... representation.
        '''
        print('convert_Keplerian_to_Alpha: Not yet implemented')
        assert False

    def convert_Alpha_to_Keplerian(self, **kwargs):
        '''
            Inverse of convert_Keplerian_to_Alpha
        '''
        print('convert_Alpha_to_Keplerian: Not yet implemented')
        assert False

    # -------------------------------------------------------------------------
    # Private Methods
    # -------------------------------------------------------------------------

    def _get_object_RADEC_from_horizons(self, objectName, obsCode, times):
        '''
            Query horizons for the RA & DEC of a *known object*
            at a sequence of times.
        '''
        print('_get_object_RADEC_from_horizons: Not yet implemented')
        assert False

    def _get_object_XYZ_from_horizons(self, objectName, obsCode, times):
        '''
            Query horizons for the BARYCENTRIC VECTOR of a *known object*
            at a sequence of times
        '''
        print('_get_object_XYZ_from_horizons: Not yet implemented')
        assert False

    def _get_orbit_RADEC(self, **kwargs):
        '''
            Propagate a specified orbit and then calculate the expected
            RA, DEC at a sequence of times.
             - wrapper around lower level method
            [[** to make life easier, demand observatory position be
              directly supplied ?? **]]

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


    def _propagate_orbit_keplerian(self, ):
        '''
            propagate an orbit using keplerian methods
             - perhaps steal import/code from cheby_checker?
             - perhaps steal import/code from neocp var-orb?

        '''
        print('_propagate_orbit_keplerian: Not yet implemented')
        assert False
        return

    def _propagate_orbit_nbody(self, ):
        '''
            propagate an orbit using nbody methods
             - perhaps steal import/code from cheby_checker?
        '''
        print('_propagate_orbit_nbody: Not yet implemented')
        assert False
        return

    def _convert_XYZ_to_RADEC(self, ):
        '''
            convert XYZ positions to apparent RA, Dec
             - perhaps steal import/code from cheby_checker?
        '''
        print('_convert_XYZ_to_RADEC: Not yet implemented')
        assert False
        return

    # -------------------------------------------------------------------------
    # Convenience funcs while developing ...
    # -------------------------------------------------------------------------

    def _interpolate_radec_for_sedna(self, times):
        '''
            '''
        # Get the look-up arrays
        JD_, RA_, Dec_ = self._radec_for_sedna()

        # Interpolate the RA & Dec at the input times
        return np.interp(times, JD_, RA_), np.interp(times, JD_, Dec_)

    def _interpolate_radec_for_101583(self, times):
        '''
            '''
        # Get the look-up arrays
        JD_, RA_, Dec_ = self._radec_for_101583()

        # Interpolate the RA & Dec at the input times
        return np.interp(times, JD_, RA_), np.interp(times, JD_, Dec_)

    def _radec_for_sedna(self,):
        '''
            '''
        return _radec_from_file(obj='sedna')

    def _radec_for_101583(self,):
        '''
            '''
        return _radec_from_file(obj='101583')


def _interpolate_radec(times, inputJRD):
    '''
        Interpolate the RA & Dec at the input times
        input:
        times - array of times for output
        inputJRD - tuple of JD_, RA_ and Dec_ of data for interpolation
    '''

    JD_, RA_, Dec_ = inputJRD

    # Interpolate the RA & Dec at the input times
    return np.interp(times, JD_, RA_), np.interp(times, JD_, Dec_)


def _radec_from_file(obj='sedna'):  # Thus also works for 101583
    JD_, RA_, Dec_ = np.genfromtxt(os.path.join(DATA_DIR, obj + '_ephem.txt'),
                                   usecols=(0, 1, 2), unpack=True)
    return JD_, RA_, Dec_


# End of file
