from __future__ import print_function

from ase.optimize.optimize import Optimizer
import numpy as np
from scipy.optimize import minimize

from ase.parallel import rank

from ase.optimize.gpmin.gp import GaussianProcess
from ase.optimize.gpmin.kernel import SquaredExponential
from ase.optimize.gpmin.prior import ConstantPrior


class GPMin(Optimizer, GaussianProcess):
    def __init__(self, atoms, restart=None, logfile='-', trajectory=None, Prior=None,
                 master=None, noise=0.005, weight=1., prior='maximum',
                 scale=0.4, force_consistent=None, batch_size=5, update=False):
        '''Beautiful doc string'''

        # DEFINING OTHER ATRIBUTES
        self.nbatch = batch_size
        self.prior = prior
        self.update_hp = update

        self.function_calls = 1
        self.force_calls = 0

        Optimizer.__init__(self, atoms, restart, logfile,
                           trajectory, master, force_consistent)

        if Prior is None:
            if prior == 'init':
                self.update_prior = False
            else:
                self.update_prior = True
            constant = self.atoms.get_potential_energy(
                force_consistent=self.force_consistent)
            Prior = ConstantPrior(constant)
        else:
            self.update_prior = False

        Kernel = SquaredExponential()
        GaussianProcess.__init__(self, Prior, Kernel)

        self.x_list = []  # Training set features
        self.y_list = []  # Training set targets

        self.set_hyperparams(np.array([weight, scale, noise]))

    def acquisition(self, r):
        e = self.predict(r)

        return e[0], e[1:]

    def update(self, r, e, f):

        # update the training set
        self.x_list.append(r)
        f = f.reshape(-1)
        y = np.append(np.array(e).reshape(-1), -f)
        self.y_list.append(y)

        if self.update_prior:
            if self.prior == 'average':
                av_e = np.mean(np.array(self.y_list)[:, 0])
                self.Prior.set_constant(av_e)
            elif self.prior == 'maximum':
                max_e = np.max(np.array(self.y_list)[:, 0])
                self.Prior.set_constant(max_e)

        # update hyperparams
        if self.update_hp and self.function_calls % self.nbatch == 0 and self.function_calls != 0:
            self.fit_to_batch()

        # build the model
        self.train(np.array(self.x_list), np.array(self.y_list))

    def relax_model(self, r0):

        result = minimize(self.acquisition, r0, method='L-BFGS-B', jac=True)

        if result.success == True:
            return result.x
        else:
            self.dump(np.array(self.x_list), np.array(self.y_list))
            raise RuntimeError(
                "The minimization of the acquisition function has not converged")

    def fit_to_batch(self):
        '''Fit hyperparameters and collect exception'''
        try:
            self.fit_hyperparameters(np.asarray(
                self.x_list), np.asarray(self.y_list))
        except Exception:
            pass

    def step(self, f):
        '''This method will be run in loop by the run method of the 
           parent class Optimizer.'''

        # Get atomic quantities
        atoms = self.atoms
        r0 = atoms.get_positions().reshape(-1)
        e0 = atoms.get_potential_energy(force_consistent=self.force_consistent)
        #f0 = f.reshape(-1)
        self.update(r0, e0, f)

        r1 = self.relax_model(r0)
        self.atoms.set_positions(r1.reshape(-1, 3))
        e1 = self.atoms.get_potential_energy(
            force_consistent=self.force_consistent)
        f1 = self.atoms.get_forces()

        self.function_calls += 1
        self.force_calls += 1

        count = 0
        while e1 >= e0:  # Or better line search condition

            self.update(r1, e1, f1)
            r1 = self.relax_model(r0)

            self.atoms.set_positions(r1.reshape(-1, 3))
            e1 = self.atoms.get_potential_energy(
                force_consistent=self.force_consistent)
            f1 = self.atoms.get_forces()

            self.function_calls += 1
            self.force_calls += 1

            if self.converged(f1):
                break

            count += 1
            if count == 30:
                raise RuntimeError('A descent model could not be built')

    def dump(self, x_train, y_train):
        '''Overwrite the method dump to be able to append data to the training set
           as it is being generated'''
        if rank == 0 and self.restart is not None:
            np.savez(self.restart, X=x_train, Y=y_train, prior=self.constant,
                     hyperparams=self.hyperparams)

    def load(self):
        '''load training set '''
        return np.load(self.restart)

    def read(self):
        data = self.load()
        self.x_list = data['X']
        self.y_list = data['Y']

