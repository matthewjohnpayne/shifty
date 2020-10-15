import MPC_library
import numpy as np 
import sys, os

class DETECTION(object):
    """
        This is intended to hold most/all of the data for a single detection
         - Perhaps not the most efficient structure, but hopefully clear
        """
    
    def __init__(self):

        # To hold data directly downloaded from horizons
        # All are for "faking" observational data ...
        #  ... or for checking the derived quantities below (in arrows, etc)
        self.horizons_cart     = None # Assume heliocentric ecliptic
        self.horizons_elem     = None # Assume heliocentric ecliptic
        self.horizons_radec    = None # Assume topocentric equatorial 
        self.horizons_unit     = None # Assume topocentric equatorial 

        # These would be required generally, ...
        # ... but I'm going to be supplying in a slightly odd way I guess 
        self.time               = None # Assume TDB (may need a conversion prior to loading) 
        self.observatory_posn   = None # Assume heliocentric ecliptic
        self.unit               = None # Assume topocentric ecliptic

        # These are for an assumed transformation distance (e.g. 2.5AU) 
        self.heliocentric_posn  = None # Assume heliocentric ecliptic
        self.heliocentric_unit  = None # Assume heliocentric ecliptic
                                       
        
        

    # ----------------------------------------------
    def ADD_HORIZONS_CARTESIAN(self, cart):
        ''' # Assume heliocentric ecliptic '''
        self.horizons_cart = cart 

    def ADD_HORIZONS_ELEMENTS(self, elem):
        ''' # Assume heliocentric ecliptic '''
        self.horizons_elem = elem
    
    def ADD_HORIZONS_RADEC(self, radec):
        ''' # Assume topocentric equatorial  '''
        self.horizons_radec = radec

    def ADD_HORIZONS_UV(self, unit):
        ''' # Assume topocentric equatorial  '''
        self.horizons_unit = unit


    # ----------------------------------------------
    def ADD_TIME(self,time):
        self.time = time

    def ADD_OBSERVATORYPOSN(self,observatory_posn):
        ''' # Assume heliocentric ecliptic '''
        self.observatory_posn = observatory_posn

    def ADD_UNIT(self,unit):
        ''' # Assume topocentric ecliptic '''
        self.unit = unit

    # ----------------------------------------------
    def ROUGH_HELIO_CONVERSION(self, r_hel ):
        ''' This does the heliocentric tranformation for the assumed radius
            Note MJP deleted any iterative stuff for LTT because I'm just doing a fixed distance
            # Assume heliocentric ecliptic
        '''
        rho_r_p, rho_r_m = self.adjust_position(r_hel, self.unit , self.observatory_posn )

        self.heliocentric_posn = np.array(rho_r_p[1])
        self.heliocentric_unit = self.heliocentric_posn/np.linalg.norm(self.heliocentric_posn)
        
        
    def adjust_position(self, r, rho_target, re):
        ''' This returns the topocentric distances and new heliocentric
        position vectors to the target, given the assumed distance
        r and the position vector of the observatory,re.
        '''
        rho_x, rho_y, rho_z = rho_target
        xe, ye, ze = re
        Robs = np.sqrt(xe * xe + ye * ye + ze * ze)
        cos_phi = -(rho_x * xe + rho_y * ye + rho_z * ze) / Robs
        phi = np.arccos(cos_phi)
        sin_phi = np.sin(phi)

        xx2 = r*r - Robs*sin_phi * Robs*sin_phi

        if xx2 < 0:
            None, None

        xx = np.sqrt(xx2)
        yy = Robs * cos_phi

        rho_p = yy + xx

        # This could be done with numpy arrays
        x_p = xe + rho_p*rho_x
        y_p = ye + rho_p*rho_y
        z_p = ze + rho_p*rho_z

        rho_m = yy - xx

        # This could be done with numpy arrays    
        x_m = xe + rho_m*rho_x
        y_m = ye + rho_m*rho_y
        z_m = ze + rho_m*rho_z

        return (rho_p, (x_p, y_p, z_p)), (rho_m, (x_m, y_m, z_m))



class TRACKLET(DETECTION):

    def __init__(self):

        # Make detection methods available
        DETECTION.__init__(self) 
        # To hold any detections
        self.dets_ = []

        # To hold data generated as an intrinsic part of the ...
        # ... transformation and clustering routine
        self.ArrowDict__ = {}

    # ----------------------------------------------
    def ADD_DETECTION(self,detection):
        self.dets_.append( detection )

    # ----------------------------------------------
    def ADD_ARROW(self,g_gdot, ref_vec, ref_t):
        self.ArrowDict__[g_gdot] = ARROW(self, g_gdot, ref_vec, ref_t)

    # ----------------------------------------------
    def GET_TIMES(self):
        self.times = np.array([det.time for det in self.dets_])

class ARROW(object):

    def __init__(self, TRACKLET , g_gdot, ref_vec, ref_t):
        # The data to work on ...
        self.T     = TRACKLET
        # The Transformation to make ...
        self.g     = g_gdot[0]
        self.gdot  = g_gdot[1]
        # The ref direction & time ...
        # ... vec is the reference direction : assume its in ecliptic coordinates 
        self.ref_vec = ref_vec/np.linalg.norm(ref_vec)
        self.ref_mat = self.xyz_to_proj_matrix(ref_vec)
        self.ref_t   = ref_t 
        # Other shit for convenience ...
        self.GM=MPC_library.Constants.GMsun
        # Do the damn transformation
        self.alphaVec = None
        self.do_transformation()

    def do_transformation(self):

        # --- N.B.: As per MJH, the first few calculations are INDEP of g/gdot -------------------------
        # Take the detection heliocentric unit vectors 
        # (from rough transform assuming 2.5AU) and rotate to projection coords 
        # Do this for each of the detections in the tracklet

        #self.theta_ = [np.dot(self.ref_mat, det.heliocentric_unit ) for det in self.T.dets_ ]
        self.theta_ = [np.dot(self.ref_mat, det.unit ) for det in self.T.dets_ ]
   
        # Rotate the observatory posns for each of the detections in the tracklet
        self.rot_observatory_posn = [np.dot(self.ref_mat, det.observatory_posn ) for det in self.T.dets_ ]

        # Light travel times
        self.LTTs_ = np.array([ vec[2] / MPC_library.Constants.speed_of_light  for vec in  self.rot_observatory_posn ])
        # ----------------------------------------------------------------------------------------------

        # Fit the tracklets ...
        alpha, alpha_dot, beta, beta_dot = self.fit_tracklet_func()

        self.alphaVec = np.array([alpha, alpha_dot, beta, beta_dot, self.g, self.gdot])

    def fit_tracklet_func(self): #, t_ref, g, gdot, v):
        # Here's a version that incorporates radial gravitational
        # acceleration
        # MJP: I believe that MJH's ``obs'' looked like: jd_tdb, dlt, theta_x, theta_y, theta_z, xe, ye, ze
        #                                                 0      1     2        3        4       5   6   7
        
        t_emit = self.T.times - self.LTTs_
        self.t_emit = t_emit - self.ref_t

        acc_z = -self.GM*self.g*self.g
        self.fac =[(1.0 + self.gdot*t + 0.5*self.g*acc_z*t*t - self.g*xyzE[2]) for xyzE, t in zip(self.rot_observatory_posn, self.t_emit)]

        A = np.vstack([self.t_emit, np.ones(len(self.t_emit))]).T
        self.A=A

        self.x = [theta[0]*f + xyzE[0]*self.g for theta, xyzE , f in zip(self.theta_, self.rot_observatory_posn, self.fac)]
        alpha_dot, alpha = np.linalg.lstsq(A, self.x)[0]

        
        self.y = [theta[1]*f + xyzE[1]*self.g for theta, xyzE , f in zip(self.theta_, self.rot_observatory_posn, self.fac)]
        beta_dot, beta = np.linalg.lstsq(A, self.y)[0]

        return (alpha, alpha_dot, beta, beta_dot)

    def xyz_to_proj_matrix(self, r_ref):
        '''This routine returns the 3-D rotation matrix for the 
        given reference vector.
        # mat is a rotation matrix that converts from ecliptic
        # vectors to the projection coordinate system.
        # The projection coordinate system has z outward,
        # x parallel to increasing ecliptic longitude, and
        # y northward, making a right-handed system
        '''
        x_ref, y_ref, z_ref = r_ref
        r = np.sqrt(x_ref*x_ref + y_ref*y_ref + z_ref*z_ref)
        lon0 = np.arctan2(y_ref, x_ref)
        lat0 = np.arcsin(z_ref/r)
        slon0 = np.sin(lon0)
        clon0 = np.cos(lon0)
        slat0 = np.sin(lat0)
        clat0 = np.cos(lat0)

        mat = np.array([[-slon0, clon0, 0], 
                        [-clon0*slat0, -slon0*slat0, clat0], 
                        [clon0*clat0, slon0*clat0, slat0 ]])

        return mat





# This rotation is taking things from equatorial to ecliptic
rot_mat = MPC_library.rotate_matrix(-MPC_library.Constants.ecl)
def equatorial_to_ecliptic(v, rot_mat=MPC_library.rotate_matrix(-MPC_library.Constants.ecl)):
    return np.dot(rot_mat, v.reshape(-1, 1)).flatten()





# Sort of psuedo-code to keep track of required steps ... 
''' 
(i) Get various input data from horizons, etc, and do what ever transformations
(ii) Put all data into DETECTIONs & do the rough-helio-conversion
(iii) Step through the nightly time-chunks and add pairs of DETECTIONS to function as TRACKLETS
(iv) For a/each g-gdot, do ADD_ARROW to each TRACKLET

'''

def convertObs80(line):
    objName   = line[0:5]
    provDesig = line[5:12]
    disAst    = line[12:13]
    note1     = line[13:14]
    note2     = line[14:15]
    dateObs   = line[15:32]
    RA        = line[32:44]
    Dec       = line[44:56]
    mag       = line[65:70]
    filt      = line[70:71]
    obsCode   = line[77:80]
    return objName, provDesig, disAst, note1, note2, dateObs, RA, Dec, mag, filt, obsCode
