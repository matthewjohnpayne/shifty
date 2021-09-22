# -*- coding: utf-8 -*-
# shifty/shifty/transformer.py

'''
Classes / methods for transforming to/from abg.
Provides methods to
-
'''

# -----------------------------------------------------------------------------
# Third party imports
# -----------------------------------------------------------------------------
import os
import sys
import getpass
# from datetime import datetime
# import copy
import numpy as np
from astropy import units as u
from astropy import constants
from astropy.time import Time
from astropy.coordinates import SkyCoord

# -----------------------------------------------------------------------------
# Any local imports
# -----------------------------------------------------------------------------
# Different machines set up differently ...
# ... adding paths to force stuff to work while developing
if getpass.getuser() in ['matthewjohnpayne']:
    sys.path.append('/Users/matthewjohnpayne/Envs/mpcvenv/')
else:
    pass

sys.path.append(os.path.dirname(os.path.dirname(
                os.path.realpath(__file__))))
# from shifty.known import Known

# -----------------------------------------------------------------------------
# Define some constants
# -----------------------------------------------------------------------------
# I'm not sure why pylint is unable to find the members insice constants and u
# but the warnings are annoying, so I'm turning those off.
# pylint: disable=no-member
speed_of_light = constants.c
speed_of_light_auPday = constants.c.to(u.au / u.day).value
au_km = u.au.to(u.km)

# -----------------------------------------------------------------------------
# Various class definitions for *data import * in shifty
# -----------------------------------------------------------------------------


class Transformer():
    '''
    (1)Keeps track of all the different coordinate systems
    (2)Converts abg to xyz
    (3)Converts abg to theta

    methods:
    --------

    main public method:
    -------------------
    abg2theta()
    get_light_travel_times()

    '''
    # Turning off some stupid syntax-checker warnings:
    # pylint: disable=too-many-instance-attributes
    # Why on earth should an object only have 7 attributes?
    # That seems dumb. Turning off this warning.
    # pylint: disable=too-few-public-methods
    # I get why an object should have at least two public methods in order to
    # not be pointless, but it's an annoying warning during dev. Turning off.

    def __init__(self, times=np.array([2459479.0, 2459480.0]),
                 time0=2459479.0, radec0=np.array([0.0, +0.0]),
                 verbose=False):
        '''
        inputs:
        -------
        time0       - float            - JD date (UTC) for reference time
        times       - array of floats  - JD date (UTC) for observations
        radec0      - list/array len=2 - Reference RA and Dec (in degrees)
        verbose     - bool             - Print extra stuff if True
        '''
        # Make sure time is an array or list
        times = np.array([times]) if (isinstance(times, int) |
                                      isinstance(times, float)) else times
        self.times = times
        self.time0 = time0
        self.verbose = verbose
        self.proj_mat = self._radec_to_proj_matrix(radec0)
        self.obs_code = None
        self.method = None

    def __call__(self, abg=np.array([0, 0, 50., 0, 0, 0]), obs_code='500@-95',
                 method='MPC'):
        '''
        input:
        abg         - array length 6  - array containing alpha, beta, gamma,
                                       alpha-dot, beta-dot, gamma-dot.
        obs_code    - string          - Observatory code.
                                        If using Horizons, note that Horizons
                                        uses some weird ones sometimes,
                                        like "500@-95" for Tess.
        method      - 'JPL' or 'MPC'  - I don't really like using MPC code
                                        as it wouldn't be publicly available(?).
                                        But having a no-internet option would
                                        also be good.
        '''
        self.obs_code = obs_code
        self.method = method
        return self.abg2theta(abg)

    def abg2theta(self, abg):
        '''
        Converts input abg to a theta vector at time dtime from reference time.
        inputs:
        -------
        abg    - array length 6 - array containing alpha, beta, gamma,
                                  alpha-dot, beta-dot, gamma-dot.
        dtime  - float          - time after reference time.
        '''
        dtime = (self.times - self.time0) / 365.24
        grav = self._grav_pert(abg)
        xyz_E = self._xyz_observer()
        x_E, y_E, z_E = xyz_E.T
        if self.verbose:
            print(xyz_E)
        num_x = abg[0] + abg[3] * dtime + abg[2] * grav[0] - abg[2] * x_E
        num_y = abg[1] + abg[4] * dtime + abg[2] * grav[1] - abg[2] * y_E
        denominator = 1 + abg[5] * dtime + abg[2] * grav[2] - abg[2] * z_E
        theta_x = num_x / denominator
        theta_y = num_y / denominator
        return np.array([theta_x, theta_y]).T

    def get_light_travel_times(self):
        '''
        Calculate the Light Travel Time.
        I'm not actually sure where I need this.
        '''
        # Light travel times (in days?)
        LTTs = np.array([vec[2] / speed_of_light_auPday
                         for vec in self._xyz_observer()])

        return LTTs

    def _grav_pert(self, abg):
        '''
        g(t), the gravitational perturbation vector, calculated from equations
        (2), (3) and (4) of B&K 2000.
        For now, just using perturbations=0, sufficient for TNOs & t<<1 yr.
        '''
        # print("grav_pert not implemented, yet. "
        #       "Enjoy zero gravity while it lasts!")
        # print(abg)
        return np.zeros(3)  # should this be a 6 or 9 vector? (derivatives)?

    def _xyz_observer(self):
        '''
        X_E(t) vector.
        Calculates the locations of the observer relative to the reference.
        For now, just make the observer not move. Only really works if t<2 days.
        '''
        helio_obs_xyz = self._get_observatory_position()
        helio_obs_xyz0 = self._get_observatory_position(reference=True)
        rel_helio_obs_xyz = helio_obs_xyz - helio_obs_xyz0
        projection_obs_xyz = self._do_transformation(rel_helio_obs_xyz)
        return projection_obs_xyz
        #return np.zeros([12,3])

    def _get_observatory_position(self, reference=False):
        '''
        Query horizons for the observatory position at a sequence of times.
        input:
        reference   - boolean - False = use self.times
                              - True  = use self.time0 reference time
        '''
        if reference:
            args = {'times': self.time0, 'obs_code': self.obs_code,
                    'verbose': self.verbose}
        else:
            args = {'times': self.times, 'obs_code': self.obs_code,
                    'verbose': self.verbose}
        # Use Horizons if explicitly requested:
        if (len(self.obs_code) != 3) | (self.method == 'JPL'):
            return get_heliocentic_ecliptic_XYZ_horizons(**args)
        # Otherwise, try using MPC tools first.
        try:
            return get_heliocentric_ecliptic_xyz_from_MPC(**args)
        except:  # If MPC tools fail, use Horizons.
            return get_heliocentic_ecliptic_XYZ_horizons(**args)

    def _do_transformation(self, observatory_posn):
        '''
        Rotate the observatory position from ecliptic cartesian to
        projection coordinates.
        '''
        # Rotate the observatory position for each time
        rotated_observatory_posn = [np.dot(self.proj_mat, obspos)
                                    for obspos in observatory_posn]
        return np.array(rotated_observatory_posn)

    def _radec_to_proj_matrix(self, radec_ref=np.array([0., 0.])):
        '''This routine returns the 3-D rotation matrix for the
        given reference ra & dec.
        # mat is a rotation matrix that converts from ecliptic
        # vectors to the projection coordinate system.
        # The projection coordinate system has z outward,
        # x parallel to increasing ecliptic longitude, and
        # y northward, making a right-handed system
        '''
        coord_ref = SkyCoord(ra=radec_ref[0] * u.deg, dec=radec_ref[1] * u.deg,
                             distance=50 * u.au)
        # I don't actually think the distance assumed changes the final result,
        # but just in case, I'm assuming 50 au above.
        x_ref, y_ref, z_ref = coord_ref.cartesian.xyz.value
        r = np.sqrt(x_ref * x_ref + y_ref * y_ref + z_ref * z_ref)
        lon0 = np.arctan2(y_ref, x_ref)
        lat0 = np.arcsin(z_ref / r)
        slon0 = np.sin(lon0)
        clon0 = np.cos(lon0)
        slat0 = np.sin(lat0)
        clat0 = np.cos(lat0)

        mat = np.array([[-slon0, clon0, 0],
                        [-clon0 * slat0, -slon0 * slat0, clat0],
                        [clon0 * clat0, slon0 * clat0, slat0]])

        return mat


# -------------------------------------------------------------------------
# These functions really don't need to be methods, and therefore aren't.
# They are more versatile (and easier to test) as functions.
# No need to over-complicate things.
# -------------------------------------------------------------------------


def get_heliocentic_ecliptic_XYZ_horizons(times, obs_code='500',
                                          verbose=False):
    '''
    Query horizons for the ecliptic heliocentric
    observatory position at a sequence of times.

    input:
    obs_code    - string
                - Note that Horizons uses some weird ones sometimes,
                  like "500@-95" for Tess.
    times       - array of JD times (UTC)
    '''
    from astroquery.jplhorizons import Horizons
    times_AP = Time(times, format='jd', scale='utc')
    # convert times to tdb
    times_tdb = times_AP.tdb.value
    horizons_query = Horizons(id='10', location=obs_code,
                              epochs=times_tdb, id_type='id')
    horizons_vector = horizons_query.vectors(refplane='ecliptic')
    helio_OBS_jpl = 0 - np.array([horizons_vector['x'], horizons_vector['y'],
                                  horizons_vector['z']]).T
    if verbose:
        print('No verbosity implemented yet, sorry')
    return helio_OBS_jpl


def get_heliocentric_equatorial_xyz_from_MPC(times, obs_code='500',
                                             verbose=False):
    '''
    Get the heliocentric EQUATORIAL vector coordinates of the
    observatory at the time jd_utc.
    '''
    # MPC_library imported here, as it is an optional dependancy
    from mpcpp import MPC_library as mpc
    obsCodes = mpc.Observatory()
    helio_OBS_equ = []
    # Make sure time is an array or list
    times_utc = np.array([times]) if (isinstance(times, int) |
                                      isinstance(times, float)) else times
    for jd_utc in times_utc:
        hom = obsCodes.getObservatoryPosition(obsCode=obs_code, jd_utc=jd_utc,
                                              old=False)
        helio_OBS_equ.append(hom)
        if verbose:
            print('MPC XYZ:')
            print(f'Heliocentric position of observatory: {hom} au\n')

    return np.array(helio_OBS_equ)


def get_heliocentric_ecliptic_xyz_from_MPC(times, obs_code='500',
                                           verbose=False):
    '''
    Get the heliocentric ECLIPTIC vector coordinates of the
    observatory at the time jd_utc.

    input:
    obs_code    - string
    times       - JD time (UTC)
    '''
    helio_OBS_equ = get_heliocentric_equatorial_xyz_from_MPC(times, obs_code,
                                                             verbose)
    helio_OBS_ecl = []
    for hequ in helio_OBS_equ:
        helio_OBS_ecl.append(equatorial_to_ecliptic(hequ))
    return np.array(helio_OBS_ecl)


def equatorial_to_ecliptic(input_xyz, ecliptic_to_equatorial=False):
    '''
    Convert an cartesian vector from mean equatorial to mean ecliptic.
    ecliptic_to_equatorial=True converts backwards, from ecliptic to equatorial.
    input:
        input_xyz              - np.array length 3
        ecliptic_to_equatorial - boolean
    output:
        output_xyz - np.array length 3
    '''
    # MPC_library imported here, as it is an optional dependancy
    from mpcpp import MPC_library as mpc
    direction = -1 if ecliptic_to_equatorial else +1
    rotation_matrix = mpc.rotate_matrix(-mpc.Constants.ecl * direction)
    output_xyz = np.dot(rotation_matrix, input_xyz.reshape(-1, 1)).flatten()
    return output_xyz


def eq2ecl(rad, decd, reverse=False):
    '''
    Function for quickly converting from equatorial (RA, Dec)
    to ecliptic latitude and longitude (l, b).
    '''
    epsilon = np.radians(-23.43928 if reverse else 23.43928)
    ra = np.radians(rad)
    dec = np.radians(decd)
    sine = np.sin(epsilon)
    cose = np.cos(epsilon)
    sina = np.sin(ra)
    cosa = np.cos(ra)
    sind = np.sin(dec)
    cosd = np.cos(dec)
    sinb = cose * sind - sine * cosd * sina
    tanl = (cose * cosd * sina + sine * sind) / (cosd * cosa)
    return np.arctan(tanl), np.arcsin(sinb)

# END
