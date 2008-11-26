import numpy as npy
from ase.transport.tools import dagger

# ## add ended by ---- if no ,oneline, if modify, comment the line before

class LeadSelfEnergy:
    #---implemented by Jingzhe Chen-------
    #nume    // number of energy points 
    #energy     // energy point in integral path,list
    #wgp     // integral weight factor in integral path,list
    #fgp     // the Fermi Distribution factor in integral path,list
    #nres    // the num of residure in the integral countour,list
    #sigma_mm     // the stored leadselfenergy cooresponding to Egp, list of matrix
    #bias          // the shift value for hamiltonian
    #kpts         //kpoints,list
        
    
    conv = 1e-8 # Convergence criteria for surface Green function
    
    
    def __init__(self, hsk_dii, hsk_dij, eta=1e-4):
        self.hk_ii, self.sk_ii = hsk_dii # onsite principal layer
        self.hk_ij, self.sk_ij = hsk_dij # coupling between principal layers
                                        # coupling to the central region
        self.nbf = self.hk_ii.shape[-1] # nbf for the scattering region
        self.eta = eta
        
        self.nkpts = self.hk_ii.shape[1]
        self.sigmak_mm = [npy.empty((self.nbf, self.nbf), complex)] * \
                                                                    self.nkpts
        self.sigma = npy.empty([self.nbf, self.nbf], complex)
        self.nume = [0] * self.nkpts
        self.nres = [0] * self.nkpts
        self.egp = []
        self.wgp = []
        self.fgp = []
        for i in range(self.nkpts):
            self.egp.append([])
            self.wgp.append([])
            self.fgp.append([])
        self.bias = 0


    def __call__(self, energy, pkpt=0):
        """Return self-energy (sigma)"""
        if energy != None:
            z = energy - self.bias 
            self.sigma = npy.empty([self.nbf, self.nbf], complex)
            tau_im = npy.empty([1, self.nbf, self.nbf], complex )
            tau_mi = npy.empty([1, self.nbf, self.nbf], complex )            
            a_im = npy.empty([1, self.nbf, self.nbf], complex )            

            tau_im = z * self.sk_ij[pkpt] - self.hk_ij[pkpt]                
            sgf = self.get_sgfinv(energy, pkpt)
            a_im = npy.linalg.solve(sgf, tau_im)
            tau_mi = z * dagger(self.sk_ij[pkpt]) - dagger(self.hk_ij[pkpt])
            self.sigma = npy.dot(tau_mi, a_im)
        return self.sigma

    def get_lambda(self, energy, pkpt=0):
        """Return the lambda (aka Gamma) defined by i(S-S^d).

        Here S is the retarded selfenergy, and d denotes the hermitian
        conjugate.
        """
        ##
        energy = npy.array(energy)
        sigma_mm = self(energy, pkpt)
        return 1.j * (sigma_mm - npy.resize(dagger(sigma_mm), sigma_mm.shape))

    def get_sgfinv(self, energy, pkpt=0):
        """The inverse of the retarded surface Green function""" 
        z = energy - self.bias 
        delta = self.conv + 1
       

        v_00 = npy.empty([self.nbf, self.nbf], complex)
        v_01 = npy.empty([self.nbf, self.nbf], complex)
        v_10 = npy.empty([self.nbf, self.nbf], complex)
        v_11 = npy.empty([self.nbf, self.nbf], complex)        
        
        v_00 = z * dagger(self.sk_ii[pkpt]) - dagger(self.hk_ii[pkpt])
        v_11 = v_00.copy()
        v_10 = z * self.sk_ij[pkpt] - self.hk_ij[pkpt]
        v_01 = z * dagger(self.sk_ij[pkpt]) - dagger(self.hk_ij[pkpt])
        n = 0
        while delta > self.conv:
            a = npy.linalg.solve(v_11, v_01)
            b = npy.linalg.solve(v_11, v_10)
            v_01_dot_b = npy.dot(v_01, b)
            v_00 -= v_01_dot_b
            v_11 -= npy.dot(v_10, a) 
            v_11 -= v_01_dot_b
            v_01 = -npy.dot(v_01, a)
            v_10 = -npy.dot(v_10, b)
            delta = npy.abs(v_01).max()
            n += 1
        return v_00  
    
    
    def set_bias(self, bias):
        self.bias = bias

        
    
