from __future__ import print_function

from ase.optimize.optimize import Optimizer
import numpy as np
from numpy import linalg as la
from scipy.optimize import minimize

from ase.parallel import rank, barrier, paropen
import pickle

from .lightGP import LightGaussianProcess
from .kernel import SquaredExponential
from .prior import ConstantPrior


class LightBayes(Optimizer, LightGaussianProcess):
   def __init__(self, atoms, restart = None, logfile = '-', trajectory = None, prior_function = None,
                      master = None, report = None, noise = 0.025, weight  = 5.,
                      scale = 0.4, force_consitent = None, nsample = 4, fit = False):

      '''Beautiful doc string'''

      #DEFINING OTHER ATRIBUTES
      self.nsample = nsample

      self.InitialPrior = True
      self.fit = fit

      self.function_calls = 1
      self.force_calls = 0
      D = len(atoms)*3

      Optimizer.__init__(self, atoms, restart, logfile, trajectory, master, force_consitent)

      if prior_function is None:
         constant = self.atoms.get_potential_energy(force_consistent = self.force_consistent)
         Prior = ConstantPrior(constant)
         prior_function = Prior.prior

      Kernel = SquaredExponential(D)
      LightGaussianProcess.__init__(self, prior_function, Kernel)


      self.x_list = [] # Training set features
      self.y_list = [] # Training set targets

      self.set_hyperparams(np.array([weight,scale, noise]))

      if report is not None:
         self.report = Reporter(report)
         self.report.prior(weight, scale, noise)
      else:
         self.report = None

   def adquisition(self, r):
      e = self.predict(r)

      return e[0]  , e[1:]

   def update(self, r, e, f):

      #update the training set
      self.x_list.append(r)
      f = f.reshape(-1)
      y = np.append(np.array(e).reshape(-1), -f)
      self.y_list.append(y)

      #build the model
      self.train(np.array(self.x_list),np.array(self.y_list))

   def relax_model(self, r0):

      result = minimize(self.adquisition, r0, method='L-BFGS-B', jac = True)

      if result.success == True:
         return result.x
      else:
         self.dump(np.array(self.x_list),np.array(self.y_list))
         raise RuntimeError("The minimization of the acquisition function has not converged")



   def step(self, f):
      '''This method will be run in loop by the run method of the
         parent class Optimizer.'''

      if self.report is not None:
         self.report.step(self.nsteps)

      #Get atomic quantities
      atoms = self.atoms
      r0 = atoms.get_positions().reshape(-1)
      e0 = atoms.get_potential_energy(force_consistent = self.force_consistent)
      #f0 = f.reshape(-1)
      self.update(r0,e0,f)

      # Fit hyperparameters in the begining
      if self.fit and self.InitialPrior and self.function_calls >= self.nsample:
         hp = self.fit_hyperparameters(np.asarray(self.x_list), np.asarray(self.y_list))
         if self.report is not None:
            arguments = (np.asarray(self.x_list), np.asarray(self.y_list))
            logP , DlogP = self.neg_log_likelihood(self.hyperparams, np.asarray(self.x_list),
                                                                     np.asarray(self.y_list))
            self.report.fit_report(hp[0], hp[1], -logP)
            self.InitialPrior = False


      r1 = self.relax_model(r0)
      self.atoms.set_positions(r1.reshape(-1,3))
      e1 = self.atoms.get_potential_energy(force_consistent = self.force_consistent)
      f1 = self.atoms.get_forces()

      if self.report is not None:
         f_max = max(np.sqrt(np.sum(f1**2, axis=1)))
         self.report.ev(self.function_calls, e1, f_max)

      self.function_calls +=1
      self.force_calls +=1


      count = 0
      while e1>=e0: #Or better line search condition

         self.update(r1,e1,f1)
         r1 = self.relax_model(r0)

         self.atoms.set_positions(r1.reshape(-1,3))
         e1 = self.atoms.get_potential_energy(force_consistent = self.force_consistent)
         f1 = self.atoms.get_forces()

         if self.report is not None:
            f_max = max(np.sqrt(np.sum(f1**2, axis=1)))
            self.report.ev(self.function_calls, e1, f_max)

         self.function_calls +=1
         self.force_calls +=1
         count +=1
         if count==30:
            raise RuntimeError('A descent model could not be built')

   def dump(self, x_train, y_train):
      '''Overwrite the method dump to be able to append data to the training set
         as it is being generated'''
      if rank == 0 and self.restart is not None:
         np.savez(self.restart, X = x_train, Y = y_train, prior = self.constant,
                  hyperparams = self.hyperparams)

   def load(self):
      '''load training set '''
      return  np.load(self.restart)

   def read(self):
      data = self.load()
      self.x_list = data['X']
      self.y_list = data['Y']



class Reporter():
   def __init__(self, filename):
      self.f = filename
      f = paropen(filename, 'w')
      f.write('PROFILE OF THE BAYESIAN OPTIMIZER\n')
      f.close()

   def prior(self,  weight, l, noise):
      f = paropen(self.f , 'a')
      f.write('Prior setups:\n ------------\n')
      f.write('Squared Exponential kernel with\n')
      f.write('\tweight: %.3f\n' % weight)
      f.write('\tscale: %.3f\n' % l)
      f.write('The Gaussian Process is regularized adding noise %.1E\n' % noise)
      f.write('\n')
      f.close()

   def step(self, nsteps):
      f = paropen(self.f, 'a')
      f.write('Step %d:\n'%nsteps)
      f.close()

   def ev(self, nfev, e, f_max):
      f = paropen(self.f, 'a')
      f.write('\tEvaluation #%d\tEnergy: %.6f\tForce %.6f\n' % (nfev,e,f_max))
      f.close()

   def fit_report(self, weight, l, likelihood):
      f = paropen(self.f, 'a')
      f.write('\n\tFIT OF THE HYPERPARAMETERS:\n')
      f.write('\tlog marginal likelihood: %.2E\n' % likelihood)
      f.write('\tweight: %.3f\tscale: %.3f\n\n' % (weight, l))
      f.close()


#--------------------------------------------------------

if __name__=="__main__":
   from ase import Atoms
   from ase.calculators.emt import EMT
   import random

   random.seed(32)

   a = 3.8 #angstrom
   e = 5 #noise

   def r(eps):
     ''' decorator'''
     return random.uniform(-eps/2.,eps/2.)

   atoms = Atoms('3Au', [(-0.3,0.,-1.1),(a+2+r(e), 0.-r(e), a+r(e)),(5.+r(e), a+r(e), a+r(e)) ])
   atoms.set_calculator(EMT())
   op = Bayes(atoms, logfile = 'bayes.log', trajectory='bayes.traj', gradient_step = 0.005, lcb = -0.005, noise = 1e-3, nsample = 0)
   op.run(0.05, 20)
   #print(op.params)
