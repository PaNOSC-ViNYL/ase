import numpy as np
from itertools import product

class Prior():
   '''Base class for all priors for the bayesian optimizer.
 
      The __init__ method and the prior method are implemented here.
      Each child class should implement its own potential method, 
      that will be called by the prior method implemented here.

      When used, the prior should be initialized outside the optimizer 
      and the prior object should be passed as a function to the optimizer.
   '''

   def __init__(self):
      '''Basic prior implementation. 
      '''
      pass

   def prior(self, x):

      '''Prior function that is called by the Gaussian Process in
         the optimizer. '''

      try:
         d = x.shape[1]
         n = x.shape[0]

         return np.hstack([self.potential(x[i,:]) for i in range(n)])

      except IndexError:
         return self.potential(x)


class ZeroPrior(Prior):
   def __init__(self):
      Prior.__init__(self)

   def potential(self, x):
      return np.zeros(x.shape[0]+1)


class ConstantPrior(Prior):
   def __init__(self, constant):
      self.constant = constant
      Prior.__init__(self)

   def potential(self, x):
      d = x.shape[0]
      output = np.zeros(d+1)
      output[0] = self.constant
      return output


class CalculatorPrior(Prior):

   def __init__(self, atoms, calculator):
      
      Prior.__init__(self)
      self.atoms = atoms.copy()
      self.atoms.set_calculator(calculator)

   def potential(self, x):
      self.atoms.set_positions(x.reshape(-1,3))
      V = self.atoms.get_potential_energy(force_consistent = True)
      gradV = -self.atoms.get_forces().reshape(-1)
      return np.append(np.array(V).reshape(-1), gradV)


