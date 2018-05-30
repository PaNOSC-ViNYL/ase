from __future__ import print_function
import numpy as np
import numpy.linalg as la
from itertools import product as iterproduct

class Kernel():
   def __init__(self):
      pass
   
   def set_params(self, params):
      pass

   def kernel(self, x1, x2):
      '''Kernel function to be fed to the Kernel matrix'''
      pass

   def K(self, X1, X2):
      '''Compute the kernel matrix '''
      return np.block([[self.kernel(x1,x2) for x2 in X2] for  x1 in X1])



class SE_kernel(Kernel):
   '''Squared exponential kernel without derivatives '''
   def __init__(self):
      Kernel.__init__(self)
      
   def set_params(self, params):
      '''Set the parameters, for example, after optimization of hyperparameters of the 
      Gaussian Process via maximizing the marginal likelihood.
      params should be a 1D array of shape (D+12, ) where D is the dimension,
      containing the bias of the Kernel as a first entry, 
      the weight as a second and the scales as the remaining'''
      
      self.weight = params[0]
      self.l = params[1]

   def SquaredDistance(self, x1, x2):
      '''Returns the norm of x1-x2 using diag(l) as metric '''

      return np.sum((x1-x2)* (x1-x2))/self.l**2

   def kernel(self,x1, x2): 
      ''' This is the squared exponential function'''
      return self.weight**2*np.exp(-0.5 * self.SquaredDistance(x1,x2))

   def dK_dweight(self, x1, x2):
      '''Derivative of the kernel respect to the weight '''
      return 2*self.weight*np.exp(-0.5 * self.SquaredDistance(x1,x2))

   def dK_dl(self, x1, x2):
      '''Derivative of the kernel respect to the scale'''
      return self.kernel*la.norm(x1-x2)**2/self.l**3



class SquaredExponential(SE_kernel):
   '''Squared exponential kernel with derivatives. 
   For the formulas, see Koistinen, Dagbjartssdittir, Asgeirsson, Vehtari, Jonsson, 
   Nudged elastic band calculations accelerated with Gaussian process regression.
   section 3.

   Atributes:
   ----------------
   weight: 	float. Multiplicative constant to the exponenetial kernel
   l : 		float. Lenght scale of the squared exponential kernel

   Relevant Methods:
   ----------------
   set_params: 		Set the parameters of the Kernel, i.e. change the atributes
   kernel_function: 	squared exponential covariance function
   kernel: 		covariance matrix between two points in the manifold. 
           		   Note the inputs are arrays of shape (D,)
   K: 			kernel matrix of two sets of data.
                           Note the inputs are arrays of shape (nsamples, D)
   gradient:            Gradient of K(X,X) with respect to the parameters of the kernel
                           i.e. the hyperparameters of the Gaussian process.
   '''


   def __init__(self, dimensionality):
      self.D = dimensionality
      SE_kernel.__init__(self)
      
   def kernel_function(self, x1, x2):
      ''' This is the squared exponential function'''
      return  self.weight**2*np.exp(-0.5 * self.SquaredDistance(x1,x2))

   def kernel_function_gradient(self, x1, x2):
      '''Gradient of kernel_function respect to the second entry.
      x1: first data point
      x2: second data point'''
          
      prefactor = (x1-x2)/self.l**2
      #return prefactor * self.kernel_function(x1,x2)
      return prefactor

   def kernel_function_hessian(self, x1, x2):
      '''Second derivatives matrix of the kernel function '''
      
      P = np.outer(x1-x2, x1-x2)/self.l**2
      prefactor = (np.identity(self.D) - P) /self.l**2

      #return prefactor * self.kernel_function(x1,x2)
      return prefactor

   def kernel(self, x1 ,x2):
      '''Squared exponential kernel including derivatives. 
      This function returns a D+1 x D+1 matrix, where D is the dimension of the manifold'''
      '''
      k =np.array([[1]]) # np.asarray(self.kernel_function(x1,x2)).reshape((1,1)) 
      j2 = self.kernel_function_gradient(x1,x2).reshape(1, -1)
      j1 = -j2.T #self.kernel_function_gradient(x2,x1).reshape(-1, 1)
      h = self.kernel_function_hessian(x1, x2)
      '''
      K = np.identity(self.D+1)
      K[0,1:] = self.kernel_function_gradient(x1,x2)
      K[1:,0] = -K[0,1:]
      #K[1:,1:] = self.kernel_function_hessian(x1, x2)
      P = np.outer(x1-x2, x1-x2)/self.l**2
      K [1:, 1:] = (K[1:, 1:]-P)/self.l**2
      #return np.block([[k,j2],[j1,h]])*self.kernel_function(x1, x2)
      return K* self.kernel_function(x1, x2)

   def kernel_matrix(self, X, nsample):
      '''This is the same method than self.K for X1=X2, but using the matrix is then symmetric'''
      #rename parameters
      D = self.D
      n = nsample

      #allocate memory
      K = np.identity((n*(D+1)), dtype = float)

      #fill upper triangular:
      for i in range(0,n):
         for j in range(i+1,n):
            k = self.kernel(X[i,:], X[j,:])
            K[i*(D+1):(i+1)*(D+1),j*(D+1):(j+1)*(D+1)] = k
            K[j*(D+1):(j+1)*(D+1),i*(D+1):(i+1)*(D+1)] = k.T
         K[i*(D+1):(i+1)*(D+1),i*(D+1):(i+1)*(D+1)] = self.kernel(X[i,:], X[i,:])
         
      return K 

   def kernel_vector(self, x, X, nsample):
      '''
      #rename parameters
      D = self.D
      n = nsample

      #allocate memory
      k = np.empty((D+1, n*D + n), dtype=float)
      for i in range(n):
         k[:, i*(D+1):(i+1)*(D+1)] = self.kernel(x, X[i,:])
      return k''' #not worth it unless n>>10
      return np.hstack([self.kernel(x, x2) for x2 in X])



   #---------Derivatives--------

   def dK_dweight(self, X):
      '''Return the derivative of K(X,X) respect to the weight '''
      return self.K(X,X)*2/self.weight

   #----Derivatives of the kernel function respect to the scale ---
   def dK_dl_k(self, x1,x2):
      '''Returns the derivative of the kernel function respect to  l '''
      #return la.norm(x1-x2)**2/self.l**3 * self.kernel_function(x1,x2)
      return np.dot((x1-x2),(x1-x2))/self.l**3
 
   def dK_dl_j(self, x1,x2):
      '''Returns the derivative of the gradient of the kernel function respect to l'''
      prefactor = -2* (1 - 0.5*self.SquaredDistance(x1,x2))/self.l
      #return self.kernel_function_gradient(x1, x2)* prefactor
      #return self.kernel_function(x1,x2)* self.kernel_function_gradient(x1, x2)* prefactor
      return self.kernel_function_gradient(x1, x2)* prefactor

   def dK_dl_h(self, x1,x2):
      '''Returns the derivative of the hessian of the kernel function respect to l'''
      I = np.identity(self.D)
      P = np.outer(x1-x2,x1-x2)/self.l**2 
      prefactor = 1-0.5*self.SquaredDistance(x1,x2)

      #return -2*self.kernel_function(x1,x2)*(prefactor*(I-P) - P)/self.l**3
      return -2*(prefactor*(I-P) - P)/self.l**3

   def dK_dl_matrix(self, x1,x2):

      k = np.asarray(self.dK_dl_k(x1,x2)).reshape((1,1))
      j2 = self.dK_dl_j(x1,x2).reshape(1, -1)
      j1 = self.dK_dl_j(x2,x1).reshape(-1, 1)
      h = self.dK_dl_h(x1, x2)


      return np.block([[k,j2],[j1,h]])*self.kernel_function(x1,x2)


   def dK_dl(self, X):
      '''Return the derivative of K(X,X) respect of l '''
      
      return np.block([[self.dK_dl_matrix(x1,x2) for x2 in X] for x1 in X])

   def gradient(self, X):
      '''Computes the gradient of matrix K given the data respect to the hyperparameters
      Note matrix K here is self.K(X,X)

      returns a 2-entries list of n(D+1) x n(D+1) matrices '''

      g = [self.dK_dweight(X), self.dK_dl(X)]

      return g




if __name__ == "__main__":
   x1 = np.array([[1,2,3]])
   x2 = np.random.rand(2,2)
   kernel = SquaredExponential(2)
   #kernel.set_params(np.array([100.,1.,6.e3,6.e3,6.e3]))
   kernel.set_params(np.array([0.5,0.3]))
   #x2 = np.array([[3.], [3.000001]])
   K = kernel.K(x2, x2)
   print(K)
   K = kernel.kernel_matrix(x2,2)
   print(K) 

   def is_pos_def(x):
    j,D = x.shape
    minors = []
    for d in range(D):
       minors.append( np.linalg.det(x[:d,:d])>0)
    return np.all(minors)   

   def is_symm(x):
      n1, n2 = x.shape
      return np.all([[x[a,b] == x[b,a] for a in range(n1)] for b in range(n2)])
   
   print(is_symm(K))

   #print(np.linalg.det(K[:3,:3]))
   print(is_pos_def(K))
   #print(np.linalg.eigvals(K))
   #from scipy.linalg import cholesky
   #L = cholesky(K)
   #G = np.array(kernel.gradient(x2))
   #print(G)
