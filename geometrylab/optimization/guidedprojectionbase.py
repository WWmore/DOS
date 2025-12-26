#!/usr/bin/env python

# -*- coding: utf-8 -*-

from __future__ import absolute_import

from __future__ import print_function

from __future__ import division

from timeit import default_timer as timer

import numpy as np

from scipy import sparse

# try: ## spsolve is faster from pypardiso, but may be incompatible with other packages
#     from pypardiso import spsolve 
# except:
from scipy.sparse.linalg import spsolve


__author__ = 'Davide Pellis'


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class GuidedProjectionBase(object):

    epsilon = 0.001

    step = 1

    iterations = 1 #changed from 10 Hui note: =each running how many optimiztions+print

    threshold = 1.0e-8

    make_residual = False

    step_control = False

    fairness_reduction = 0

    verbose = True

    _weights = {}

    _N = 0

    _X = None

    _X0 = None

    _H = None

    _r = None

    _H0 = None

    _r0 = None

    _K = None

    _s = None

    _K0 = None

    _s0 = None

    _errors = {}

    _values = {}

    _messages = []

    _report = None

    _counter = 0

    _initialize = True

    _reinitialize = True

    _constant_constraints = []

    _iterative_constraints = []

    _constraints = []

    _residual = None

    __timer = 0

    __t0 = 0

    __counter = 0

    #--------------------------------------------------------------------------
    #
    #--------------------------------------------------------------------------

    @property
    def N(self):
        return self._N

    @N.setter
    def N(self, N):
        if N != self._N:
            self.reinitialize = True
        self._N = N

    @property
    def X(self):
        return self._X

    @X.setter
    def X(self, X):
        self._X = X
        self._X0 = np.copy(X)

    @property
    def max_weight(self):
        return max(list(self._weights.values()))

    @property
    def weights(self):
        return self._weights

    @property
    def initialize(self):
        return self._initialize

    @initialize.setter
    def initialize(self, bool):
        self._initialize = bool

    @property
    def reinitialize(self):
        return self._reinitialize

    @reinitialize.setter
    def reinitialize(self, bool):
        self._reinitialize = bool

    @property
    def iteration(self):
        return self._counter

    #--------------------------------------------------------------------------
    #                                 Workfolw
    #--------------------------------------------------------------------------

    def _initialization(self):
        self._reset_counter()
        self.on_initialize()

    def initialization(self):
        self._initialization()
        self.reinitialization()

    def _reinitialization(self):
        self.on_reinitialize()
        self._X = None
        self._X0 = None
        self._residual = None
        self._errors = {}
        self._counter = 0
        self.initialize_unknowns_vector()
        #self._make_errors()# Hui note: comment to unshow the residual table
        #self.make_values()
        if self.verbose:
            print('  *** Initialized ***')
            #print(self._report) # Hui note: comment to unshow the residual table

    def reinitialization(self):
        self._set_weights()
        self._set_dimensions()
        self._reinitialization()

    def _initialize_iteration(self):
        if self.initialize:
            self._initialization()
            self._initialize = False
            self._reinitialize = True
        self._set_weights()
        self._set_dimensions()
        self.reinitialize_check()
        if self.reinitialize:
            self._reinitialization()
            self._reinitialize = False
        self._constraints = []
        self._constant_constraints = []
        self.__t0 = timer()

    def _pre_iteration_update(self):
        self.balance_weights()
        self._constraints = self._constant_constraints + []

    def _post_iteration_update(self):
        self.post_iteration_update()
        self._counter += 1
        self.__counter += 1
        self.make_values()

    def _finalize(self):
        self.__timer += timer() - self.__t0
        #self._make_X0() #Huinote:should comment, otherwise X0 equal X
        #self._make_errors()# Hui note: comment to unshow the residual table
        # if self.verbose:# Hui note: comment to unshow the residual table
        #     print(self._report)
        # self.on_finalize()

    def _reset_counter(self):
        self._counter = 0
        self.__counter = 0
        self.__timer = 0

    def _make_X0(self):
        X = np.copy(self.X)
        self._X = None
        self.initialize_unknowns_vector()
        self._X = X

    #--------------------------------------------------------------------------
    #                               Initialization
    #--------------------------------------------------------------------------

    def _set_weights(self):
        self.set_weights()

    def _set_dimensions(self):
        self.set_dimensions()

    def _make_errors(self):
        self._errors = {}
        self._messages = []
        self.make_messages()
        self.make_errors()
        self._make_report()

    #--------------------------------------------------------------------------
    #                                 Overwrite
    #--------------------------------------------------------------------------

    def reinitialize_check(self):
        '''Launch here the function self.reinitialize()'''
        pass

    def on_reinitialize(self):
        pass

    def on_initialize(self):
        pass

    def on_finalize(self):
        pass

    def initialize_unknowns_vector(self):
        pass

    def set_weights(self):
        pass

    def set_dimensions(self):
        pass

    def balance_weights(self):
        pass

    def make_errors(self):
        pass

    def make_messages(self):
        pass

    def make_values(self):
        pass

    def post_iteration_update(self):
        '''Update here objects from the unknowns'''
        pass

    #--------------------------------------------------------------------------
    #                                  Add
    #--------------------------------------------------------------------------

    def add_error(self, name, mean_error, max_error, weight=None):
        self._errors[name] = (mean_error, max_error, weight)

    def add_value(self, name, value):
        self._values[name] = value

    def add_message(self, message):
        self._messages.append(message)

    def add_weight(self, name, default=1):
        self._weights[name] = default

    def add_weights(self, weights):
        self._weights = {**self._weights, **weights}

    #--------------------------------------------------------------------------
    #                                  Get
    #--------------------------------------------------------------------------

    def get_value(self, key):
        try:
            return self._values[key]
        except:
            return None

    def get_error(self, key):
        try:
            return self._errors[key]
        except:
            return None

    def get_weight(self, key):
        return self._weights[key]

    #--------------------------------------------------------------------------
    #                                  Set
    #--------------------------------------------------------------------------

    def set_weight(self, name, value):
        if name in self._weights:
            self._weights[name] = value
        else:
            out = ('{} weight does not exists!').format(name)
            raise ValueError(out)

    #--------------------------------------------------------------------------
    #                                Format
    #--------------------------------------------------------------------------

    def error_string(self, name):
        error = self.get_error(name)
        if error is None:
            return '_'
        string = '{:.4E} | {:.4E}'
        out = (string).format(error[0], error[1])
        return out

    #--------------------------------------------------------------------------
    #                                 Report
    #--------------------------------------------------------------------------

    def _out_error(self, name):
        error = self.get_error(name)
        if error is None:
            error = ('-', '-', '-')
            string = '.{:>22}: {:>11} | {:>10} | {:>10} .\n'
        else:
            string = '.{:>22}: {:>11.4E} | {:>10.4E} | {:>10.4E} .\n'
        out = (string).format(name, *error)
        return out

    def _format_message(self, message):
        out = '. ' + str(message) + '\n'
        return out

    def _split_string(self):
        n = 64
        out = '.'*n + '\n'
        return out

    def _open_string(self):
        n = 64
        out = '-'*n + '\n'
        name = 'iteration ' + str(self.__counter)
        out += '.' + ' '*((n-len(name))//2) + name + '\n'
        out += '-'*n + '\n'
        name = 'Relative Error [mean | max | weight]'
        out += '.' + ' '*((n-len(name))//2) + name + '\n'
        return out

    def _optimization_infos(self):
        out = ('. cumulative time = {:.3f}\n').format(self.__timer)
        out += ('. number of variables = {}\n').format(self.N)
        constraints = 0
        for c in self._constraints:
            constraints += c[1]
        out += ('. number of constraints = {}\n').format(constraints)
        return out

    def _make_report(self):
        out =  self._open_string()
        out += self._split_string()
        for error in self._errors:
            out += self._out_error(error)
        out += self._split_string()
        out += self._optimization_infos()
        out += self._split_string()
        for message in self._messages:
            out += self._format_message(message)
        out += self._split_string()
        self._report = out
        if self.make_residual:
            self.__make_residuals_report()
        return out

    def _make_residuals_report(self):
        if self._residual is None:
            return
        out = '. Residuals\n'
        out += self._split_string()
        O = 0
        for constraint in self._constraints:
            res = self._residual[O:O + constraint[1]]
            res = np.linalg.norm(res)
            out += ('. {} = {:.8f}\n').format(constraint[0], res)
            O += constraint[1]
        out += self._split_string()
        self._report += out

    # -------------------------------------------------------------------------
    #                                 Save
    # -------------------------------------------------------------------------

    def save_report(self, file_name):
        try:
            name = ('{}_report.txt').format(file_name)
            txt = open(name, 'w')
            txt.write(self._report)
            txt.close()
        except:
            print('Report not available!')

    # -------------------------------------------------------------------------
    #                                Build
    # -------------------------------------------------------------------------

    def add_iterative_constraint(self, H, r, name='constraint'):
        self._H = sparse.vstack((self._H, H))
        self._r = np.hstack((self._r, r))
        self._constraints.append((name, H.shape[0]))

    def add_constant_constraint(self, H, r, name='constraint'):
        self._H0 = sparse.vstack((self._H0, H))
        self._r0 = np.hstack((self._r0, r))
        self._constant_constraints.append((name, H.shape[0]))

    def add_iterative_fairness(self, K, s, name='fairness'):
        self._K = sparse.vstack((self._K, K))
        self._s = np.hstack((self._s, s))

    def add_constant_fairness(self, K, s, name='fairness'):
        self._K0 = sparse.vstack((self._K0, K))
        self._s0 = np.hstack((self._s0, s))

    def build_iterative_constraints(self):
        pass

    def build_constant_constraints(self):
        pass

    def build_iterative_fairness(self):
        pass

    def build_constant_fairness(self):
        pass

    def build_regularizer(self):
        self._R = self.epsilon ** 2 * sparse.identity(self.N)

    # -------------------------------------------------------------------------
    #                                Build
    # -------------------------------------------------------------------------

    def _build_constant_matrices(self):
        null =  np.zeros([0])
        self._H0 = sparse.coo_matrix((null,(null,null)), shape=(0,self.N))
        self._K0 = sparse.coo_matrix((null,(null,null)), shape=(0,self.N))
        self._r0 = np.array([])
        self._s0 = np.array([])
        self.build_constant_constraints()
        self.build_constant_fairness()

    def _build_iterative_matrices(self):
        null =  np.zeros([0])
        self._H = sparse.coo_matrix((null,(null,null)), shape=(0,self.N))
        self._K = sparse.coo_matrix((null,(null,null)), shape=(0,self.N))
        self._r = np.array([])
        self._s = np.array([])
        self.build_iterative_constraints()
        self.build_iterative_fairness()
        self.build_regularizer()

    def _make_residuals(self):
        H = sparse.vstack([self._H0, self._H])
        r = sparse.hstack([self._r0, self._r]).T
        r = r.toarray()
        X = sparse.csc_matrix([self.X]).transpose()
        self._residual = H.dot(X) - r

    #--------------------------------------------------------------------------
    #                                Solver
    #--------------------------------------------------------------------------

    def optimize(self):
        self._initialize_iteration()
        self._build_constant_matrices()
        X0 = np.array(self.X)
        #diff = np.linalg.norm(X0) / X0.shape[0]
        iteration = 0
        #while diff > self.threshold and iteration < self.iterations:
        for i in range(self.iterations):
            self._pre_iteration_update()
            self._build_iterative_matrices()
            H = sparse.vstack([self._H0, self._H])
            r = np.array([np.hstack((self._r0, self._r))]).T
            K = sparse.vstack([self._K0, self._K])
            s = np.array([np.hstack((self._s0, self._s))]).T
            K = 1.0/(10**(iteration * self.fairness_reduction)) * K
            s = 1.0/(10**(iteration * self.fairness_reduction)) * s
            X = sparse.csc_matrix([self.X]).transpose()
            
            #Hui changed to below since an error: spmatrix.dot has not dot
            A = H.T @ H + K.T@K + self._R
            a = H.T @ r + K.T@s + self.epsilon**2 * X
            # A = (sparse.spmatrix.dot(H.transpose(), H)
            #      + sparse.spmatrix.dot(K.transpose(), K) + self._R)
            #a = sparse.spmatrix.dot(H.transpose(), r) \
                       #+ sparse.spmatrix.dot(K.transpose(), s) \
                       #+ self.epsilon**2 * X
                         
            try:
                a = a.toarray()
            except AttributeError:
                pass
            if self.step_control and self._residual is None:
                res0 = H.dot(X) - r.toarray()
                res0 = np.linalg.norm(res0)
                self._residual = res0
                #print(res0)
            
            #try: ##Huinote: original
            X = spsolve(A,a)
            # except: ##Huinote: added
            #     from pypardiso import PyPardisoSolver
            #     solver = PyPardisoSolver()
            #     solver.iparm[1] = 1  # 设置为 1 可以减少内存使用
            #     solver.factorize(A)
            #     X = solver.solve(A, a)
            
            iteration += 1
            Xi = np.array(X)
            if self.step != 1:
                Xi = self.step*Xi + (1-self.step)*X0
            elif self.step_control:
                lam = 1
                stop = False
                r = r.toarray()
                while not stop:
                    X = lam*Xi + (1-lam)*X0
                    X = sparse.csc_matrix([X]).transpose()
                    res = H.dot(X) - r
                    res = np.linalg.norm(res)
                    if res < self._residual or lam < 1e-2:
                        self._residual = res
                        stop = True
                        Xi = lam*Xi + (1-lam)*X0
                    lam = 0.5*lam
            #diff = np.linalg.norm(X0-Xi) / X0.shape[0]
            X0 = np.array(Xi)
            self._X = np.array(Xi)
            self._post_iteration_update()
        self._finalize()
        return X0