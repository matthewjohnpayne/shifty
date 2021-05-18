import numpy as np
import sys, os
from functools import lru_cache
import MPC_library
import scipy
from collections import defaultdict
import DES_IO_Utilities as DESIO

# Useful utility functions I should make
#  (i) Load/populate a single exposure and run calc for single set of (ref_vec, ref_t, gamma, gammadot, alphadot, betadot)
# (ii) Load/populate a single exposure and run calc as for (i) but for ARRAY OF alphadot & betadot
#(iii) Load/populate multiple exposures and run calc for single set of (ref_vec, ref_t, gamma, gammadot, alphadot, betadot)
# (iv) Load/populate multiple exposures and run calc as for (iii) but for ARRAY OF alphadot & betadot

'''
    def make_and_run_XXX( exposureDataDict , searchParamsDict , detectionDataDict, DEBUG=False):

    # Check that the input data is of the right type ...
    if not (isinstance(exposureDataDict, dict) and isinstance(searchParamsDict, dict) and isinstance(detectionDataDict, dict)):
        sys.exit("\n\tDemand all of these be dictionaries: types=(%r,%r ,%r)" % (type(exposureDataDict), type(searchParamsDict), type(detectionDataDict)))
    # Check that the input dicts have the right keys in them ...
    dict_key_pairs = [
        [exposureDataDict,  ['expTime','observatory_posn']],
        [searchParamsDict,  ['ref_vec','ref_t', 'gamma','gammadot','alphadot','betadot']],
        [detectionDataDict, ['RA','DEC']]
    ]
    for dict_key_pair in dict_key_pairs:
        Dict,RequiredKeys = dict_key_pair
        if not np.all( [ key in Dict for key in RequiredKeys ]):
            sys.exit("\n\tExpected the following keys %r , but only found %r " % (RequiredKeys, list(Dict.keys()) ) )
    

    # Create the exposure object
    EXP = EXPOSURE()
    # Add the exposure-specific info
    EXP.ADD_TIME(exposureDataDict['expTime'])
    EXP.ADD_OBSERVATORYPOSN(exposureDataDict['observatory_posn'])
    # Add in the observations
    EXP.ADD_RADEC(detectionDataDict['RA'],detectionDataDict['DEC'])

    # Run the search
    EXP.transform_detections_to_dots(searchParamsDict['ref_vec'],searchParamsDict['ref_t'],
                 searchParamsDict['gamma'],searchParamsDict['gammadot'],
                 searchParamsDict['alphadot'],searchParamsDict['betadot'],
                 DEBUG=DEBUG)

    return EXP
'''

class EXPOSURE(object):
    """ This is intended to 
        (i) Hold the important exposure-level quantities
        (ii) Hold the detections
        (iii) Perform transformations using specified parameters
    """
    
    # ------------------------------------------------------------------------
    def __init__(self):
        
        # Exposure-specific quantities, ...
        self.expTime              = None # Exposure Mid Time (JD) ?
        self.observatory_posn     = None # Heliocentric coordinates of Observatory. Coord-Sys must match Detections [Equatorial]
    
        # Detection positions, untransformed & transformed...
        # ...will want these to be VECTOR quantities
        self.RA_      = None # Input RA  [Degrees], Coords = [Topocentric Equatorial]
        self.DEC_     = None # Input Dec [Degrees], Coords = [Topocentric Equatorial]
        self.unit_    = None # Unit-Vector from (RA,Dec), Coords = [Topocentric Equatorial (because calc from RA,Dec) ]
        self.thetaX_  = None # Detection posns transformed to plan perpendicular to reference vector [X-component, radians]
        self.thetaY_  = None # Detection posns transformed to plan perpendicular to reference vector [Y-component, radians]
        
        # Dict to hold DOTS
        self.DotDict = {}
        
        # Other shit for convenience ...
        self.GM=MPC_library.Constants.GMsun
    # ------------------------------------------------------------------------


    # - Populate the basic exposure quantities -------------------------------
    # --- Not all of the functions are used ... some just-in-case ------------
    def ADD_ID(self,ID):
        self.ID = ID

    def ADD_TIME(self,expTime):
        self.expTime = expTime
    
    #def ADD_OBSERVATORYPOSN(self,observatory_posn):
    #    self.observatory_posn = observatory_posn
    
    def getObservatoryPosition(self,obscode):
        # The routine below calculates the HelioEquatorialObservatoryPosn
        # in equatorial cartesian coordinates.
        # Note that if comparing to Lackners code, he may be using Barycentric
        self.observatory_posn = MPC_library.Observatories.getObservatoryPosition('%s'%obscode, self.expTime)
        # Now we convert to ECLIPTIC coordinates
        self.observatory_posn = self.equatorial_to_ecliptic(self.observatory_posn)

    def ADD_RADEC(self,RA_, DEC_):
        if not (isinstance(RA_, np.ndarray) and isinstance(DEC_, np.ndarray)):
            sys.exit("Input RA_ & DEC_ need to be np.arrays : type(RA_) = %r ,  type(DEC_) = %r " % ( type(RA_), type(DEC_)) )
        self.RA_ = RA_
        self.DEC_ = DEC_
        # We'll need the RA,DEC in vector form
        self.calculate_unit(RA_, DEC_)
    
    def ADD_UNIT(self,unit_):
        if not (isinstance(unit_, np.ndarray) and unit_.shape[-1]==3):
            sys.exit("Input unit_ needs to be np.array & needs to be of shape (X,3) : type(unit_) = %r ,  unit_.shape = %r " % ( type(unit_), unit_.shape) )
        self.unit_ = unit_

    def calculate_unit(self, RA_, DEC_):
        # This translates (RA,DEC) [assumed in degrees & assumed in Equatorial coords] into a unit vector
        x_ = np.cos(DEC_*np.pi/180.) * np.cos(RA_*np.pi/180.)
        y_ = np.cos(DEC_*np.pi/180.) * np.sin(RA_*np.pi/180.)
        z_ = np.sin(DEC_*np.pi/180.)
        # Now we convert to ECLIPTIC coordinates
        v = self.equatorial_to_ecliptic(np.array([x_,y_,z_]).T )
        # Store ECLIPTIC unit vector
        self.ADD_UNIT( v )
    
    def equatorial_to_ecliptic(self, v, rot_mat=MPC_library.rotate_matrix(-MPC_library.Constants.ecl)):
        # Convert from equatorial coordinates to ECLIPTIC coordinates
        return np.matmul(v, rot_mat.T)
    # ------------------------------------------------------------------------



    # - Functions that are NOT dependent on gamma and/or gammadot ------------
    # --- These are calculating rotns & LTTs
    # --- Applied to all dets in an exposure (for a given ref-vec & ref-time)
    def call_reference_calcs(self, ref_vecX, ref_vecY, ref_vecZ, ref_t):
        """  Wrapper around the ~5 funcs below that are only dep on ...
            ... ref_vec and/or ref_t
            
            N.B. input ref_vec is in component-form to allow cache to work
        """
        ref_vec = np.array([ref_vecX, ref_vecY, ref_vecZ])
        self.ref_vec = ref_vec
        # Calculate the rotn matrix
        rot_mat = self.calculate_rotation_matrix(ref_vec)
        
        # Rotate the observatory posns
        rot_observatory_posn = self.rotate_observatory_posn( rot_mat )
        
        # Rotate the observations to theta-coords
        thetaXYZ_ = self.rotate_the_observations_into_theta_coords(rot_mat)

        # Make them available & also return them (not sure at present which is most useful ) ...
        self.rot_mat , self.rot_observatory_posn , self.ref_t, self.thetaX_, self.thetaY_ =\
            rot_mat, rot_observatory_posn , ref_t, thetaXYZ_[:,0]/thetaXYZ_[:,2], thetaXYZ_[:,1]/thetaXYZ_[:,2]
        
        return rot_mat, rot_observatory_posn , ref_t, thetaXYZ_[:,0]/thetaXYZ_[:,2], thetaXYZ_[:,1]/thetaXYZ_[:,2]
    
    def calculate_rotation_matrix(self, ref_vec):
        ''' This routine returns the 3-D rotation matrix for the
            given reference vector.
            # rot_mat is a rotation matrix that converts from ecliptic
            # vectors to the projection coordinate system.
            # The projection coordinate system has z outward,
            # x parallel to increasing ecliptic longitude, and
            # y northward, making a right-handed system
        '''
        x_ref, y_ref, z_ref = ref_vec
        r = np.sqrt(x_ref*x_ref + y_ref*y_ref + z_ref*z_ref)
        lon0 = np.arctan2(y_ref, x_ref)
        lat0 = np.arcsin(z_ref/r)
        slon0 = np.sin(lon0)
        clon0 = np.cos(lon0)
        slat0 = np.sin(lat0)
        clat0 = np.cos(lat0)
    
        rot_mat = np.array([[-slon0, clon0, 0],
                        [-clon0*slat0, -slon0*slat0, clat0],
                        [clon0*clat0, slon0*clat0, slat0 ]])
        return rot_mat
    
    def rotate_observatory_posn(self, ref_mat):
        # Rotate the observatory posn for the exposure
        #print("rotated ob posn = ", np.dot(ref_mat, self.observatory_posn ))
        return np.dot(ref_mat, self.observatory_posn )
    

    def rotate_the_observations_into_theta_coords(self, ref_mat):
        # Rotate the observational unit-vector into projection-plane coordinates (theta)
        #print("thetaXYZ_ = ", ref_mat.dot(self.unit_.T).T)
        return ref_mat.dot(self.unit_.T).T

    # ------------------------------------------------------------------------



    # - Functions that are only dependent on gamma and/or gammadot -----------

    @lru_cache(maxsize=256)
    def calc_Alpha0_Beta0(self, gamma, gammadot):
        
        # Calculate the (relative) time of emission of the light
        self.tEmit = self.expTime - (1.0/gamma - self.rot_observatory_posn[2] )/MPC_library.Constants.speed_of_light - self.ref_t   #
        print("expTime, (1.0/gamma - self.rot_observatory_posn[2] )/C, retT, tEmit", self.expTime , (1.0/gamma - self.rot_observatory_posn[2] )/MPC_library.Constants.speed_of_light , self.ref_t, self.tEmit)
        
        # Calculate the 'scaling factors' required for the calculation of alpha0,beta0
        factor0 = 1.0 -  0.5* self.GM * gamma**3 *self.tEmit**2    - gamma*self.rot_observatory_posn[2]    + gammadot*self.tEmit
        self.factor3 = 1.0 - 0.5*self.GM * gamma**3  *self.tEmit*self.tEmit

        # Calculate & return two arrays: alpha0, beta0
        return (factor0*self.thetaX_  +  gamma*self.rot_observatory_posn[0])/self.factor3 , (factor0*self.thetaY_  +  gamma*self.rot_observatory_posn[1])/self.factor3







    # - Functions that are dep. on alpha-dot & beta-dot --------------------------------------

    # Calculate Alpha & Beta by applying the linear shifts (funcs of alphadot, betadot) ------
    def calculate_alpha_beta_SCALAR(self,alphadot, betadot):
        """ I am assuming that the necessary preceeding steps have been evaluated ...
            ... such that the following variables exist: self.Cx_ , self.Cy_, & self.tEmit
            
            This is using a SCALAR value for self.tEmit
             - Seems to be most useful internally
        """
        #alpha_, beta_ = self.Cx_ - alphadot*self.tEmit_  ,    self.Cy_ - betadot*self.tEmit_
        #self.abDict[self.abkey(alphadot, betadot)]=np.vstack((alpha_, beta_)).T
        #return np.vstack((alpha_, beta_)).T
        pass
    
    def calculate_alpha_beta_VECTOR(self,alphadot, betadot):
        """ I am assuming that the necessary preceeding steps have been evaluated ...
            ... such that the following variables exist: self.Cx_ , self.Cy_, & self.tEmit_
            
            This is using a VECTOR form for self.tEmit_
             - Seems to be the most use externally when I am jamming lots of exposures together
        """
        #alpha_, beta_ = self.Cx_ - alphadot*self.tEmit_  ,    self.Cy_ - betadot*self.tEmit_
        alpha_, beta_ = self.a0_ - alphadot*self.tau_  ,    self.b0_ - betadot*self.tau_
        
        self.abDict[self.abkey(alphadot, betadot)]=np.vstack((alpha_, beta_)).T
        return np.vstack((alpha_, beta_)).T

    def calculate_alpha_beta_VECTOR_SIMPLE(self,alphadot, betadot):
        #alpha_, beta_ = self.Cx_ - alphadot*self.tEmit_  ,    self.Cy_ - betadot*self.tEmit_
        alpha_, beta_ = self.a0_ - alphadot*self.tau_  ,    self.b0_ - betadot*self.tau_
        return np.vstack((alpha_, beta_)).T

    def abkey(self, alphadot, betadot):
        return ("%.6e" % alphadot, "%.6e" % betadot)
    # ---------------------------------------------------------------------------------------
    
    
    
    # - Functions that do CLUSTERING --------------------------------------------------------
    
    # Do Clustering at a given alpha, beta --------------------------------------------------
    def find_cluster_at_single_alphadot_betadot(self, alphadot, betadot, VECTORFUNC=True, rad=(np.pi/180.)*(15./3600.), minClusterSize=6):
        # Use "calculate_alpha_beta", then use results to build a tree
        tree = scipy.spatial.cKDTree(self.calculate_alpha_beta_VECTOR_SIMPLE(alphadot, betadot))
        
        # Query the tree
        matches = tree.query_ball_tree(tree, rad)

        # Minimum cluster size (is there a more efficient way?)
        return set( [tuple(m) for m in matches if len(m) >=minClusterSize ] )
    
    def find_clusters_over_multiple_radii_at_single_alphadot_betadot(self, alphadot, betadot, radiusArray, minClusterSize=6):
        # Use "calculate_alpha_beta", then use results to build a tree
        #print("Building tree...")
        tree = scipy.spatial.cKDTree(self.calculate_alpha_beta_VECTOR_SIMPLE(alphadot, betadot))
        # Do searches over each radius
        setList = []
        for rad in radiusArray:
            matches   = tree.query_ball_tree(tree, rad)
            setList.append( set( [tuple(m) for m in matches if len(m) >=minClusterSize ] ) )
        #print()
        return setList
    
    # Do Clustering over a *RANGE* of alphadot & betadot
    # N.B. I am also building in the possibility of multiple clusterradius searches
    def find_cluster_over_multiple_alphadot_betadot(self, alphadotArray, betadotArray, radiusArray, minClusterSize=6):
        cluster_setList_Dict = {}
        for alphadot in alphadotArray:
            for betadot in betadotArray:
                # A "setList" is returned by "find_clusters_over_multiple_radii_at_single_alphadot_betadot" ...
                # .. I put these into a dictionary ...
                cluster_setList_Dict[self.abkey(alphadot,betadot)] = self.find_clusters_over_multiple_radii_at_single_alphadot_betadot(alphadot, betadot, radiusArray, minClusterSize=minClusterSize)
        return cluster_setList_Dict

    # Max dist between any two points in arc-sec (requires that calculate_alpha_beta have been run & that self.abDict be populated)
    def calculate_max_distances(self):
        self.abDistDict = {k:np.max(scipy.spatial.distance.pdist(alpha_beta_points))*(180/np.pi)*(3600) for k,alpha_beta_points in self.abDict.items() }


    # ---------------------------------------------------------------------------------------






class DOT_SET(EXPOSURE):
    """ This is intended to be used for processing (clustering) multiple sets of EXPOSURE objects
        I know that we are going to want to process many detections from many exposures in the same manner
        So this is intended to facilitate that by performing the required calculations in a sensible manner
        
        Previously we said          'ARROWS == Transformed(TRACKLET)'
        Analagously I am now saying 'DOTS   == Transformed(DETECTION)'
        
    """

    # *** SEEMS LIKE I WILL WANT TO SLIGHTLY TWEAK THE FUNCTIONS AVAILABLE IN THE EXPOSURE CLASS ***
    # Want to ...
    # ... load exposures, perhaps into one long array
    # ... do I want to allow gamma, gamma-dot calculations at this stage, or only within exposures ?
    # ... want to be able to calculate alpha-dot, beta-dot shifts over the whole set



    def __init__(self):
        # - Allow ourselves to use exposure methods on collected sets  --------
        super().__init__()
        # - We will store results of alpha-dot,beta-dot transformations in a dictionary
        self.abDict = defaultdict(np.array)


    # - Populate the DOT_SET with EXPOSURES -------------------------------
    def ADD_EXPOSURES_calcC(self,EXPOSURES_,gamma, gammadot):
        """ For simplicity I am assuming that the EXPOSURES have had "call_reference_calcs" evaluated
            Hence the EXPOSURES should contain ... self.rot_mat , self.rot_observatory_posn , self.tEmit, self.thetaX_, self.thetaY_
        """
        #self.Cx_    = np.array([])
        #self.Cy_    = np.array([])
        self.tEmit_ = np.array([])
        self.a0_    = np.array([])
        self.b0_    = np.array([])
        self.tau_   = np.array([])
        for EXP in EXPOSURES_:
            
            # Calculate the Cx_ & Cy_ array for each of the exposures ...
            # ... and load into a common array within this DOT_SET object
            #Cx_,Cy_  = EXP.calc_C(gamma, gammadot)
            #self.Cx_ = np.append(self.Cx_ , Cx_)
            #self.Cy_ = np.append(self.Cy_ , Cy_)
            a0_,b0_   = EXP.calc_Alpha0_Beta0(gamma, gammadot)
            self.a0_ = np.append(self.a0_ , a0_)
            self.b0_ = np.append(self.b0_ , b0_)

            # It will also be convenient to have the  tEmit variable in array format ...
            #self.tEmit_ = np.append( self.tEmit_ , np.full(Cx_.shape[0], EXP.tEmit) )
            #self.tEmit_ = np.append( self.tEmit_ , np.full(a0_.shape[0], EXP.tEmit) )
            self.tau_   = np.append( self.tau_ , np.full(a0_.shape[0], EXP.tEmit/EXP.factor3) )










class EFFICIENT_EXPOSURE_SET(EXPOSURE):

    def __init__(self, name):
        # - Allow ourselves to use exposure methods
        super().__init__()
        # - Unique name : will get used for any file-saves
        self.name = name
        # Other shit for convenience ...
        self.GM=MPC_library.Constants.GMsun


    # In MJP_single_detections_A, I used a function "perform_initial_steps" to do a lot of the initial steps.
    # I want to put that functionality into an class object
    # And I want to make it as efficient as possible (both in memory & computation)
    def generate_base_quantities(self,df, SAVE=False):
    

        # (i) Get the raw exposure data
        raArray_,decArray_,expTimeArray_, observatory_posnArray_ = DESIO.get_ARRAY_data_for_many_exposures(df)
        for item in (raArray_,decArray_,expTimeArray_, observatory_posnArray_): print(len(item))
        # (ii) Calculate observational unit vectors
        # This translates (RA,DEC) [assumed in degrees & assumed in Equatorial coords] into a unit vector
        x_ = np.cos(decArray_*np.pi/180.) * np.cos(raArray_*np.pi/180.)
        y_ = np.cos(decArray_*np.pi/180.) * np.sin(raArray_*np.pi/180.)
        z_ = np.sin(decArray_*np.pi/180.)
        # Now we convert to ECLIPTIC coordinates
        rot_mat = MPC_library.rotate_matrix(-MPC_library.Constants.ecl)
        np.matmul(np.array([x_,y_,z_]).T, rot_mat.T)
        ObsUnitVectorsArray_ = self.equatorial_to_ecliptic(np.array([x_,y_,z_]).T )

        # Save ...
        if SAVE:
            self.save_basic_obs(ObsUnitVectorsArray_, observatory_posnArray_ , expTimeArray_)
        
        return ObsUnitVectorsArray_, observatory_posnArray_ , expTimeArray_



    # - Functions that are NOT dependent on gamma and/or gammadot ------------
    # --- These are calculating rotns & LTTs
    # --- Applied to all dets in an exposure (for a given ref-vec & ref-time)
    
    def calculate_theta_etc_from_arrays(self,ref_vec, ref_time, ObsUnitVectorsArray_, observatory_posnArray_ , expTimeArray_ ,  SAVE=False):
        # Calculate the rotn matrix: using EXPOSURE functionality
        rot_mat = self.calculate_rotation_matrix(ref_vec)
        
        # Rotate the observatory posns into projection-plane coordinates
        rot_observatory_posnArray_ = np.dot(rot_mat, observatory_posnArray_.T).T
        
        # Rotate the observational unit-vector into projection-plane coordinates
        projectionXYZ_ = rot_mat.dot( ObsUnitVectorsArray_.T).T
        
        # Make them available & also return them (not sure at present which is most useful ) ...
        thetaXArray_, thetaYArray_ = projectionXYZ_[:,0]/projectionXYZ_[:,2], projectionXYZ_[:,1]/projectionXYZ_[:,2]
        
        # Save ...
        if SAVE:
            self.save_theta_etc(ref_vec, ref_time, thetaXArray_, thetaYArray_, expTimeArray_, rot_observatory_posnArray_)

        return thetaXArray_, thetaYArray_, expTimeArray_, rot_observatory_posnArray_



    # - Functions that are only dependent on gamma and/or gammadot -----------
    
    def calculate_alpha0_beta0_tau_from_arrays(self,ref_vec, ref_time, gamma, gammadot, thetaXArray_, thetaYArray_, expTimeArray_, rot_observatory_posnArray_ , SAVE=False):
        # Calculate (see "calc_Alpha0_Beta0" in EXPOSURE) ...
        
        # The (relative) time of emission of the light
        emissionTimeArray_ = expTimeArray_ - (1.0/(MPC_library.Constants.speed_of_light*gamma)) - ref_time
        emissionTimeArray2_= np.square(emissionTimeArray_)

        # The 'scaling factors' required for the calculation of alpha0,beta0
        factor0 = 1.0 - 0.5*self.GM * gamma**3 *emissionTimeArray2_    - gamma*rot_observatory_posnArray_[:,2]    + gammadot*emissionTimeArray_
        factor3 = 1.0 - 0.5*self.GM * gamma**3 *emissionTimeArray2_
        
        # Calculate & return two arrays: alpha0, beta0
        alpha0Array_ = (factor0*thetaXArray_  +  gamma*rot_observatory_posnArray_[:,0])/factor3
        beta0Array_  = (factor0*thetaYArray_  +  gamma*rot_observatory_posnArray_[:,1])/factor3
        tauArray_    = emissionTimeArray_ / factor3

        # Save ...
        if SAVE:
            self.save_alpha0_beta0_tau(ref_vec, ref_time, gamma, gammadot, alpha0Array_, beta0Array_, tauArray_)

        return alpha0Array_, beta0Array_, tauArray_



    # - Functions that are dep. on alpha-dot & beta-dot --------------------------------------
    def calculate_alpha_beta_from_arrays(self,ref_vec, ref_time, gamma, gammadot, alphadot, betadot, alpha0Array_, beta0Array_, tauArray_,  SAVE=False):
        # Calculate
        alphaArray_, betaArray_ = alpha0Array_ - alphadot*tauArray_  ,    beta0Array_ - betadot*tauArray_
        # Save ...
        if SAVE:
            self.save_alpha_beta(ref_vec, ref_time, gamma, gammadot, alphadot, betadot, alphaArray_, betaArray_)
        return alphaArray_, betaArray_


    # - Functions for clustering ------------------------------------------------------------
    def find_cluster_at_single_alphadot_betadot_from_arrays(self,ref_vec, ref_time, gamma, gammadot, alphadot, betadot, alpha0Array_, beta0Array_, tauArray_,  SAVE=False, rad=(np.pi/180.)*(15./3600.), minClusterSize=6):
        # Use "calculate_alpha_beta_from_arrays", then use results to build a tree
        tree = scipy.spatial.cKDTree( np.vstack((self.calculate_alpha_beta_from_arrays(ref_vec,ref_time, gamma, gammadot, alphadot, betadot, alpha0Array_, beta0Array_, tauArray_,SAVE=False))).T )
        
        # Query the tree
        matches = tree.query_ball_tree(tree, rad)
        
        # Get only those above a minimum cluster size (is there a more efficient way?)
        matchSet = set( [tuple(m) for m in matches if len(m) >=minClusterSize ] )
        # Save ...
        if SAVE:
            self.save_cluster(ref_vec, ref_time, gamma, gammadot, alphadot, betadot, "SINGLE", matchSet)
        
        #return matches
        return (matches, matchSet)


    def abkey(self, alphadot, betadot):
        return ("%.6e" % alphadot, "%.6e" % betadot)
    def gabkey(self, gamma,alphadot, betadot):
        return ("%.6e" % gamma, "%.6e" % alphadot, "%.6e" % betadot)
    
    
    
    
    # Grid-based "clustering"
    def grid_search_over_multiple_alphadot_betadot(self,ref_vec, ref_time, gamma, gammadot, alphadotArray, betadotArray, alpha0Array_, beta0Array_, tauArray_,  SAVE=False, width=2*np.pi/(180*3600), minClusterSize=5):
        D = {}
        for alphadot in alphadotArray:
            for betadot in betadotArray:
                D[self.abkey(alphadot, betadot)] = self.grid_search(ref_vec, ref_time, gamma, gammadot, alphadot, betadot, alpha0Array_, beta0Array_, tauArray_,  SAVE=False, width=width, minClusterSize=minClusterSize)
        return D
                
    
    def grid_search(self,ref_vec, ref_time, gamma, gammadot, alphadot, betadot, alpha0Array_, beta0Array_, tauArray_,  SAVE=False, width=2*np.pi/(180*3600), minClusterSize=5):
        alphaArray_, betaArray_ = self.calculate_alpha_beta_from_arrays(ref_vec,ref_time, gamma, gammadot, alphadot, betadot, alpha0Array_, beta0Array_, tauArray_,SAVE=False)
        # Make into integers
        X = ((alphaArray_-alphaArray_.min()) / width).astype(int)
        Y = ((betaArray_-betaArray_.min()) / width).astype(int)
        # Make a unique long integer
        N=1000000
        XY = N*X+Y
        # Find the unique items, and count them
        unique_categories, unique_cat_Counts = np.unique(XY, return_counts=True)
        # This finds the counts that are >= some critical number of counts
        ind1 = np.where(unique_cat_Counts >= minClusterSize)[0]
        big_unique_cats = unique_categories[ind1]

        # This grabs/generates the 9 neighboring "keys"
        neighbors = self.OneToMany(big_unique_cats)

        # Check which of the "neighbors" is actually populated
        AltKeysPresent = unique_categories[np.isin(unique_categories, neighbors)]

        # Individual objects associated with "big_unique_cats"
        # Data looks like ...
        #   356046200578 (array([2282, 2308, 2314, 2322, 2327]),)
        #   345689191124 (array([5445, 5639, 5644, 5712, 5934]),)
        #   347266208893 (array([3970, 3977, 3987, 4332, 4357, 4364]),)
        indD = { b:np.nonzero(np.isin(XY,n[np.isin(n,AltKeysPresent)]))     for b,n in zip(big_unique_cats, neighbors) }
        return indD
    
    
    def OneToMany(self,keys_):
        A = np.around(keys_,-6)
        B = keys_-A
        A = A/1000000
        return np.array([ (1000000)*(A+i)+B+j for j in range(-1,2,1) for i in range(-1,2,1) ]).astype(int).T

    
    
    
    
    
    
    
    
    
    
    
    # - Functions that are for loading/saving ------------------------------------------------
    def save_cluster(self,ref_vec, ref_time, gamma, gammadot, alphadot, betadot, TYPE, matchSet):
        filename = self.filenameformatting([ref_vec, ref_time, gamma, gammadot, alphadot, betadot, TYPE])
        np.savez(filename, matchSet=matchSet)
    def get_cluster_from_file(self,ref_vec, ref_time, gamma, gammadot, alphadot, betadot, TYPE):
        filename = self.filenameformatting([ref_vec, ref_time, gamma, gammadot, alphadot, betadot])
        file     = np.load(filename)
        return file['matchSet']

    def save_alpha_beta(self,ref_vec, ref_time, gamma, gammadot, alphadot, betadot, alphaArray_, betaArray_):
        filename = self.filenameformatting([ref_vec, ref_time, gamma, gammadot, alphadot, betadot])
        np.savez(filename, alpha=alphaArray_, beta=betaArray_)
    def get_alpha_beta_from_file(self,ref_vec, ref_time, gamma, gammadot, alphadot, betadot):
        filename = self.filenameformatting([ref_vec, ref_time, gamma, gammadot, alphadot, betadot])
        file     = np.load(filename)
        return file['alpha'], file['beta']
    
    def save_alpha0_beta0_tau(self,ref_vec, ref_time, gamma, gammadot, alpha0Array_, beta0Array_, tauArray_):
        filename = self.filenameformatting([ref_vec, ref_time, gamma, gammadot])
        np.savez(filename, alpha0=alpha0Array_, beta0=beta0Array_, tau=tauArray_)
    def get_alpha0_beta0_tau_from_file(self,ref_vec, ref_time, gamma, gammadot):
        filename = self.filenameformatting([ref_vec, ref_time, gamma, gammadot])
        file     = np.load(filename)
        return file['alpha0'], file['beta0'], file['tau']
    
    def save_theta_etc(self,ref_vec, ref_time,  thetaXArray_, thetaYArray_, expTimeArray_, rot_observatory_posnArray_):
        filename = self.filenameformatting([ref_vec, ref_time])
        np.savez(filename, thetaX=thetaXArray_, thetaY=thetaYArray_, expTime=expTimeArray_, rot_observatory_posn=rot_observatory_posnArray_)
    def get_theta_etc_from_file(self,ref_vec, ref_time):
        filename = self.filenameformatting([ref_vec, ref_time])
        file     = np.load(filename)
        return file['thetaX'], file['thetaY'], file['expTime'], file['rot_observatory_posn']
    
    def save_basic_obs(self,ObsUnitVectorsArray_, observatory_posnArray_ , expTimeArray_):
        filename = self.filenameformatting([])
        np.savez(filename, ObsUnitVectors=ObsUnitVectorsArray_, observatory_posn=observatory_posnArray_, expTime=expTimeArray_)
    def get_basic_obs_from_file(self):
        filename = self.filenameformatting([])
        file     = np.load(filename)
        return file['ObsUnitVectors'], file['observatory_posn'], file['expTime']

    def filenameformatting(self,paramList, STEM='MJPSAVES/'):
        names = ['ref_vec', 'ref_time', 'gamma', 'gammadot', 'alphadot', 'betadot', 'TYPE']
        str = STEM + self.name + "__"
        if len(paramList) > 0 :
            D__ = {names[n]:p for n,p in enumerate(paramList)}
            str = str +  "v_{v0:.6f}_{v1:.6f}_{v2:.6f}__t_{t:.6f}__".format(v0 = D__['ref_vec'][0] , v1 = D__['ref_vec'][1] , v2 = D__['ref_vec'][2] , t = D__['ref_time'] )
        if len(paramList) > 2:
            str = str + "g_{g:.6e}__gd_{gd:.6e}__".format(g = D__['gamma'] , gd = D__['gammadot'] )
        if len(paramList) > 4:
            str = str + "ad_{ad:.6e}__bd_{bd:.6e}__".format(ad = D__['alphadot'] , bd = D__['betadot'])
        if len(paramList) > 6:
            str = str + TYPE + "__"
        print(str)
        return str + ".npz"


