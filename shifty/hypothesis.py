class OrbitHypothesis():
    '''
        Container for pixel-shifts and alpha-dot, beta-dot, gamma, gamma-dot
    '''

    def __init__(self, shift_x, shift_y , alpha_dot, beta_dot, gamma, gamma_dot):
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
        (d_min=100, d_max=300) or (alpha, beta, gammas), 
        ImageDataSet
        
        Methods:
        --------
        get_alpha_beta_gamma()
        get_hypotheses() => list of `OrbitHypothesis`

    '''
    
    def __init__(self,
                 ImageDataSet,
                 gamma_list,
                 gamma_dot_list = [0.] ):
        
        self.gamma_list = gamma_list
    

    def _get_default_alphadot_betadot_ranges(self, gamma):
        '''
            Set the default search range as a func of gamma
            [[This will be from the max bound velocity ]]
        '''
        alpha_dot_min, alpha_dot_max, beta_dot_min, beta_dot_max  = 0,1,0,1
        return alpha_dot_min, alpha_dot_max, beta_dot_min, beta_dot_max

    def generate_hypotheses():
        '''
        Will enumerate all possible alpha-dot, beta-dot, gamma (and gamma-dot?) values that will be investigated
        For each individual alpha-dot, beta-dot, gamma, it will generate the appropriate list of "shifts"
        '''

        hypotheses = []
        
        for gamma in self.gamma_list:
            alpha_dot_list, beta_dot_list = self._get_default_alphadot_betadot_ranges(gamma)
        
            for alpha_dot in alpha_dot_list:
                for beta_dot in beta_dot_list:
    
                    # Do something to generate shifts
                    # - will presumably use ImageDataSet
                    
                    # Pass shifts to OrbitHypothesis object ...
                    hypotheses.append( OrbitHypothesis(shift_x, shift_y , alpha_dot, beta_dot, gamma, gamma_dot) )


        return hypotheses
