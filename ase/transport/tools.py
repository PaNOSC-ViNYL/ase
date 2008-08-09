import numpy as npy 
la = npy.linalg
from math import sqrt

def tri2full(M,UL='L'):
    """UP='L' => fill upper triangle from lower triangle
       such that M=M^d"""
    nbf = len(M)
    if UL=='L':
        for i in range(nbf-1):
            M[i,i:] = M[i:,i].conjugate()
    elif UL=='U':
        for i in range(nbf-1):
            M[i:,i] = M[i,i:].conjugate()


def dagger(matrix, copy=1):
    # First change the axis: (Does not allocate a new array)
    matrix_t = matrix.swapaxes(0,1)
    if copy: # Allocate space for new array
        matrix_t = matrix_t.conj()    
    else:    # The array of matrix is used for output
        matrix_t.imag *= -1
        #npy.multiply(matrix_conj.imag,-1,matrix_conj.imag)
    return matrix_t

def lambda_from_self_energy(self_energy,copy=1):
    if copy:
        lambda_lead = self_energy.copy()
    else:
        lambda_lead = self_energy
    self_energy_t = self_energy.swapaxes(0,1) #no new array
    lambda_lead.real -= npy.array(self_energy_t.real) #no new array
    lambda_lead.imag += npy.array(self_energy_t.imag) #no new array
    lambda_lead *= 1.0j  #no new array

    return lambda_lead

def rotate_matrix(h,u):
    return npy.dot(u.T.conj(),npy.dot(h,u))

def get_subspace(a_ii,index_j):
    """Get the subspace spanned by the basis function listed in index_j"""
    return npy.take(npy.take(a_ii,index_j,axis=0),index_j,axis=1)

def normalize_rot(c,s):
    """normalize column vectors so that <c[:,i]|s|c[:,i]> = 1 """
    for i in xrange(c.shape[0]):
        v = c[:,i]
        norm = 1.0/sqrt(npy.dot(v.conj(),npy.dot(s,v)))
        c[:,i] *= norm

def subdiagonalize(h_ii,s_ii,index_j):
    nb = h_ii.shape[0]
    nb_sub = len(index_j)
    h_sub_jj = get_subspace(h_ii,index_j)
    s_sub_jj = get_subspace(s_ii,index_j)
    e_j,v_jj = la.eig(la.solve(s_sub_jj,h_sub_jj))
    normalize_rot(v_jj,s_sub_jj) #normalize: <v_j|s|v_j> = 1
    permute_list = npy.argsort(e_j.real)
    e_j = npy.take(e_j,permute_list)
    v_jj = npy.take(v_jj,permute_list,axis=1) # 
    #setup transformation matrix
    c_ii = npy.identity(nb,npy.complex)
    for i in xrange(nb_sub):
        for j in xrange(nb_sub):
            c_ii[index_j[i],index_j[j]] = v_jj[i,j]

    h1_ii = rotate_matrix(h_ii,c_ii)
    s1_ii = rotate_matrix(s_ii,c_ii)

    return h1_ii,s1_ii,c_ii,e_j

def cutcoupling(h,s,index_n):
    for i in index_n:
        s[:,i] = 0.0
        s[i,:] = 0.0
        s[i,i] = 1.0
        Ei = h[i,i]
        h[:,i] = 0.0
        h[i,:] = 0.0
        h[i,i] = Ei


#def lambda_from_self_energy2(self_energy):
#    return 1.0j * (self_energy - dagger(self_energy))

#def LambdaFromSelfEnergy(selfenergy,copy_selfenergy=1):
#        import copy
#        # Calculates lambda = i[sigma-sigma^dagger]
#        if copy_selfenergy: # Allocate space for a new array
#                lambda_lead=copy.copy(selfenergy)
#        else:
#                lambda_lead=selfenergy
#                
#        # Using that lambda = i [ (sigma.real-sigma.real^T)+
#        #                        i(sigma.imag+sigma.imag^T)]    
#        selfenergy_t=Numeric.swapaxes(selfenergy,0,1)
        # Numeric.array array performs a copy of the uncontiguous array
#        Numeric.subtract(lambda_lead.real,Numeric.array(selfenergy_t.real),lambda_lead.real)
#        Numeric.add(lambda_lead.imag,Numeric.array(selfenergy_t.imag),lambda_lead.imag)
#        Numeric.multiply(lambda_lead,complex(0,1),lambda_lead)
#        return lambda_lead

#def FermiDistribution(energy,kBT=0.0):
#        if kBT==0.0:
#            return 0.5*(1.0-Numeric.sign(energy))
#        value,sign=abs(energy/kBT),Numeric.sign(energy)
#        exponent=sign*min(value,1.0e2)
#        return 1.0/(1.0+Numeric.exp(exponent))
#
#def WriteToNetCDFFile(filename,x,y):
#        from Scientific.IO.NetCDF import NetCDFFile
#        file=NetCDFFile(filename,'w')
#        file.createDimension('Npoints',len(x))
#        myvara=file.createVariable('xvalues',Numeric.Float,('Npoints',))
#        myvara[:]=x
#        myvarb=file.createVariable('yvalues',Numeric.Float,('Npoints',))
#        myvarb[:]=y
#        file.sync()
#        file.close()

#def ReadFromNetCDFFile(filename):
#        from Scientific.IO.NetCDF import NetCDFFile
#        from Dacapo import NetCDF
#        file=NetCDFFile(filename,'r')
#        x=NetCDF.Entry(name='xvalues').ReadFromNetCDFFile(file).GetValue()
#        y=NetCDF.Entry(name='yvalues').ReadFromNetCDFFile(file).GetValue()
#        return x,y

#def RotateArray(a,U):
#        from ASE.Utilities.Wannier.HamiltonianTools import RotateMatrix
#        b=Numeric.zeros(a.shape,Numeric.Complex)
#        for i in range(len(a)):
#            b[i]=RotateMatrix(a[i],U)
#        return b
#
#
#def Convolution(f,g,zero_index,axis=0):
#        # Calculates the convolution (f * g)(x) = (1/(2*pi))\int dy f(y)g(x-y).
#        # This is the Fourier transform of f(t)g(t).
#        # If f and g are matrices the integral is performed along the specified axis.
#        # 'zero_index' is the index number of the zero point.
#        # 'time_translate' translates the two functions relative to each other by 
#        # the specified number og grid points.
#        import FFT
#        f_fft=FFT.fft(f,axis=axis)
#        g_fft=FFT.fft(g,axis=axis)
#        product_fft=f_fft*g_fft
#        pre=1.0/(2*Numeric.pi)
#        return pre*TranslateAlongAxis0(FFT.inverse_fft(product_fft,axis=0),-zero_index)
#
#
#def DirectConvolution(f,g,grid,de):
#        # Calculates (1/pi)\int dy f(y)g(x-y) directly (no FFT).
#        # 'grid' is in integer points. 'de' is energy spacing. 
#        s=[len(grid)]+list(f.shape[1:])
#        h=Numeric.zeros(s,Numeric.Complex)
#        for i in range(len(grid)):
#                f_translate=TranslateAlongAxis0(f,-grid[i])
#                h[i]=Numeric.sum(f_translate*g)*de 
#        return h/Numeric.pi
#
#def AntiConvolution(f,g,zero_index,axis=0):
#        # Calculates the anti-convolution (f ** g)(x) = (1/(2*pi))\int dy f(x+y)g(y).
#        # This is the Fourier transform of f(t)g(-t).
#        # If f and g are matrices the integral is performed along the specified axis.
#        # 'zero_index' is the index of the zero point.
#        import FFT
#        f_fft=FFT.fft(f,axis=axis)
#        g_fft=Numeric.conjugate(FFT.fft(Numeric.conjugate(g),axis=axis))
#        product_fft=f_fft*g_fft
#        pre=1.0/(2*Numeric.pi)
#        return pre*TranslateAlongAxis0(FFT.inverse_fft(product_fft,axis=0),zero_index)
#
#def HilbertTransform(f,energy):
#	from ASE.Transport.Hilbert import Hilbert
#	return Hilbert(f, nfft=len(energy), kerneltype='Simple')
#
#        # Below is old and obsolete code
#
#        # f is in energy space
#        f_fft=FFT.fft(f,axis=0) 
#        f_fft[0]=Numeric.zeros(f_fft[0].shape,Numeric.Complex)
#        dim=f[0].shape
#        
#        one_over_w_time=Numeric.reshape(Numeric.sign(energy)*(-1.0j),[len(energy),1,1])
#        t1=Numeric.repeat(one_over_w_time,repeats=dim[0],axis=1)
#        t2=Numeric.repeat(t1,repeats=dim[1],axis=2)
#########################################        return FFT.inverse_fft(f_fft*t2,axis=0)
