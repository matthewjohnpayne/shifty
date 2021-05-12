'''
Methods to propagate known orbits for shifty.py
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
    Class to deal with propagatation of known orbits

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

    For development only:
    _interpolate_radec_for_sedna()
    _interpolate_radec_for_101583()
    _radec_for_sedna()
    _radec_for_101583()
    '''

    def __init__(self, **kwargs):
        super().__init__()
        if 'obs_code' in kwargs and 'times' in kwargs and \
           (('object_name' in kwargs) or ('orbit' in kwargs)):
            self.get_known_RADEC(**kwargs)

    # -------------------------------------------------------------------------
    # Public Methods
    # -------------------------------------------------------------------------
    def get_known_RADEC(self, **kwargs):
        '''
        General method to get RADEC for specified object/orbit
        - calls specific sub-method as necessary:
           _get_object_RADEC_from_horizons if "object_name" supplied,
           _get_orbit_RADEC if "orbit" is supplied [not functional yet])
        '''
        if 'obs_code' not in kwargs or 'times' not in kwargs:
            sys.exit('get_known_RADEC() always requires at least '
                     '*obs_code* and an array of *times* to be specified')
        else:
            self.times = kwargs['times']
            self.obs_code = kwargs['obs_code']

        if 'object_name' in kwargs:
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

    def _get_object_RADEC_from_horizons(self, object_name, obs_code, times,
                                        object_type='smallbody'):
        '''
        Query horizons for the RA & DEC of a
        *known object* at a sequence of times.
        input:
        object_name - string
        obs_code    - string
                    - Note that Horizons uses some weird ones sometimes,
                      like "500@-95" for Tess.
        times       - array of times (JD)
        object_type - string
                    - Usually the default "smallbody" is fine, 
                      but for some objects, like natsats, it's neccessary.

        '''
        horizons_query = Horizons(id=object_name, location=obs_code,
                                  epochs=times, id_type=object_type)
        horizons_ephem = horizons_query.ephemerides(extra_precision=True)
        self.RA, self.Dec = np.array([horizons_ephem['RA'],
                                      horizons_ephem['DEC']])

    def _get_object_XYZ_from_horizons(self, object_name, times,
                                      object_type='smallbody',
                                      plane='earth'):
        '''
        Query horizons for the BARYCENTRIC VECTOR of a
        *known object* at a sequence of times.
        input:
        object_name - string
        times       - array
        object_type - string
                    - Usually the default "smallbody" is fine, 
                      but for some objects, like natsats, it's neccessary.
        plane       - string
                    - 'earth' = equatorial
                    - 'ecliptic' = ecliptic
        '''
        horizons_query = Horizons(id=object_name, location='500@0',
                                  epochs=times, id_type=object_type)
        horizons_vector = horizons_query.vectors(refplane=plane)
        self.XYZ = np.array([horizons_vector['x'], horizons_vector['y'],
                              horizons_vector['z']])

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

    def _interpolate_radec_for_sedna(self, times, obs_code='C57'):
        '''
        Interpolate the RA & Dec at the input times.
        '''
        return _interpolate_radec(times,
                                  _radec_from_file(obj='Sedna',
                                                   obs_code=obs_code))

    def _interpolate_radec_for_101583(self, times, obs_code='C57'):
        '''
        Interpolate the RA & Dec at the input times.
        '''
        return _interpolate_radec(times,
                                  _radec_from_file(obj='Sedna',
                                                   obs_code=obs_code))

    def _radec_for_sedna(self, obs_code='C57'):
        '''
        Get RA & Dec at hourly intervals from file.
        '''
        return _radec_from_file(obj='Sedna', obs_code=obs_code)

    def _radec_for_101583(self, obs_code='C57'):
        '''
            '''
        return _radec_from_file(obj='101583', obs_code=obs_code)


def _interpolate_radec(times, inputJRD):
    '''
    Interpolate the RA & Dec at the input times.
    input:
    times - array of times for output
    inputJRD - tuple of JD_, RA_ and Dec_ of data for interpolation
    '''

    JD_, RA_, Dec_ = inputJRD

    # Interpolate the RA & Dec at the input times
    return np.interp(times, JD_, RA_), np.interp(times, JD_, Dec_)


def _radec_from_file(obj='Sedna', obs_code='C57'):
    '''
    Read JD, RA & Dec from file for a given object and obs_code.
    input:
    obj      - string - object name
    obs_code - string - observatory code. 
    '''
    if obs_code=='500@-95':
        obs_code = 'C57'
    filename = obj + '_ephem_' + obs_code + '.txt'
    JD_, RA_, Dec_ = np.genfromtxt(os.path.join(DATA_DIR, filename),
                                   delimiter=(17, 5, 13, 13),
                                   usecols=(0, 2, 3), unpack=True)
    return JD_, RA_, Dec_


def _xyz_interp(times, inputJRD):
    '''
    Interpolate the XYZ at the input times
    input:
    times - array of times for output
    inputJXYZ - tuple of JD_, X_, Y_ and Z_ of data for interpolation
    '''

    JD_, X_, Y_, Z_ = inputJRD

    # Interpolate the RA & Dec at the input times
    return (np.interp(times, JD_, X_), np.interp(times, JD_, Y_),
            np.interp(times, JD_, Z_))


def _xyz_from_file(obj='Sedna'):
    '''
    Read JD, X, Y and Z from file for a given object and obs_code.
    input:
    obj      - string - object name
    '''
    filename = obj + '_vector.txt'
    JD_, X_, Y_, Z_ = np.genfromtxt(os.path.join(DATA_DIR, filename),
                                    delimiter=(17, 66, 23, 3, 23, 3, 23),
                                    usecols=(0, 2, 4, 6), unpack=True)
    return JD_, X_, Y_, Z_


# End of file
