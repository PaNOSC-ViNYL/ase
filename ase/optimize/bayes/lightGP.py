from __future__ import print_function
from kernel import SquaredExponential

import numpy as np
import numpy.linalg as la
from scipy.optimize import minimize
from scipy.linalg import solve_triangular
from scipy.linalg import cho_factor, cho_solve

from prior import ZeroPrior, ConstantPrior



class LightGaussianProcess():
   def __init__(self, prior = None, kernel=None):
      
      if kernel == None:
         self.kernel = SquaredExponential()
      else:
         self.kernel = kernel

      if prior == None:
         Prior = ZeroPrior()
         self.prior = Prior.prior
      else:
         self.prior = prior


   def set_hyperparams(self, params):
      self.hyperparams = params
      self.kernel.set_params(params[:-1])
      self.noise = params[-1]

   def train(self, X, Y, noise = None):
      '''Given a set of observations, X, Y, gradY, compute the K matrix
      of the Kernel given the data and its cholesky factorization such that
      this needs only to run once if new data is added

      inputs: 
      X: observations(i.e. positions). numpy array with shape: nsamples x D
      Y: targets (i.e. energy). numpy array with shape (nsamples, D+1) '''
 
      if noise is not None:
        self.noise = noise 	#Set noise atribute to a different value      

      self.X = X.copy()		 #Store the data in an atribute
      #self.y = Y.flatten()       #Reshape the observations
      K = self.kernel.kernel_matrix(X,len(X))     #Compute the kernel matrix
      K[range(K.shape[0]), range(K.shape[0])] += self.noise**2
       
      self.m = self.prior(X)

      self.L, self.lower = cho_factor(K, lower = True, check_finite=True)
      self.a = Y.flatten() - self.m 
      cho_solve((self.L, self.lower),self.a, overwrite_b=True, check_finite=True) 

   def predict(self, x):
      '''Given a trained and fitted gaussian process, it predicts the value and the 
      derivative of the target at x
      It returns f and V:
      f : prediction: [y, grady]
      V : Covariance matrix. Its diagonal is the variance of each component of f 


      TO BE FINISHED: explain what gradient is'''
      
      #x = x.reshape(1,-1)
      n = self.X.shape[0]
      k = self.kernel.kernel_vector(x, self.X, n)
      
      f = self.prior(x) + np.matmul(k,self.a)         

      return f


   def neg_log_likelihood(self, p , *args):
      '''Negative logarithm of the marginal likelihood.
         Nice documentation should be written here '''

      X, Y = args		
      self.kernel.set_params(np.append(p, self.noise))
      self.train(X, Y)

      y = Y.flatten()
 
      #Compute log likelihood
      logP = -0.5* np.dot(y-self.m , self.a) -np.sum(np.log(np.diag(self.L )))-X.shape[0]*0.5*np.log(2*np.pi)  

      #Gradient of the loglikelihood
      #grad = np.array(self.kernel.gradient(X))
      grad = self.kernel.gradient(X)      

      #iL = la.inv(self.L)
      #iK = np.matmul(iL.T,iL)
      a = self.a.reshape(1,-1)
      #vectorizing the derivative of the log likelyhood
      #D_P_input =np.array( [np.matmul( np.matmul(self.a.T, self.a), g) for g in grad])
      D_P_input =np.array( [np.matmul(a.T, np.matmul(a, g)) for g in grad])
      #D_P_input = np.matmul( np.outer(self.a, self.a), grad)
      D_complexity = np.array([cho_solve((self.L, self.lower) ,  g) for g in grad])
         
      #DlogP=np.array([0.5*np.trace( np.matmul( np.outer(self.a,self.a) -iK, g)) for g in grad])
      DlogP = 0.5 * np.trace(D_P_input - D_complexity, axis1=1, axis2=2)
      return -logP , -DlogP 

   def fit_hyperparameters(self, X, Y):

      '''Given a set of observations, X, Y, gradY; optimize the hyperparameters 
      of the Gaussian Process maximizing the marginal log-likelihood.
      This method calls TRAIN, so the atributes X, L and a are filled during the 
      process, there is no need to call the TRAIN method again.
      The method also sets the parameters of the Kernel to their optimal value at
      the end of execution

      inputs: 
      X: observations(i.e. positions). numpy array with shape: nsamples x D
      Y: targets (i.e. energy). numpy array with shape (nsamples, )
      gradY: derivative of target ( i.e. forces).
             numpy array with shape (nsamples,D) '''
     
      params = np.copy(self.hyperparams)[:-1]
      arguments = (X, Y)
      result = minimize(self.neg_log_likelihood, params, args = arguments, 
                         method = 'L-BFGS-B', jac = True)

      if result.success==False:
          print(result)
          raise NameError("The Gaussian Process could not be fitted.")
      else: # result.sucess==True :)
         self.hyperparams = np.append(result.x.copy(), self.noise)
         #print(self.hyperparams)
      self.set_hyperparams(self.hyperparams)
      return self.hyperparams



if __name__ == "__main__":
   
   GP = GaussianProcess()
   GP.set_hyperparams(np.array([ 1.0, 2.0, 1e-3]))
   X = np.array([[0.1,0.2,300],[0.4,0.1,1.], [300, 5., 8.]])
   Y = np.array([1000, 2.000,-1000])
   gradY = np.array([[0.,0.,0.],[1., 1., 1.],[-15,0.,-1]])
   observations = []
   for i in range(len(Y)):
      observations.append(np.block([Y[i], gradY[i]]).reshape(-1))
   observations = np.asarray(observations)
   print(observations)
   
   GP.train(X,observations, 1e-6)
   f, V = GP.predict(np.array([3,1,2]))
   print('Predicted energy', f[0])
   print('Predicted force', f[1:])
   print('Predicted variance', np.diag(V))
   
   logP, dlogP = GP.neg_log_likelihood(np.array([1., 2.]), X, observations)
   print(logP, dlogP)
   p = GP.fit_hyperparameters(X,observations)
   print(p) 
