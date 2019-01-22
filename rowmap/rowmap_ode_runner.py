#!/usr/bin/python
# -*- coding: utf-8 -*-
# Author: Andreas Buttenschoen
from __future__ import print_function

from scipy.integrate._ode import IntegratorBase
import numpy as np
import rowmap._rowmap as _rowmap

import re
import warnings

""" SciPy integrator wrapper for RowMap implemented in rowmap.f90 """
class rowmap(IntegratorBase):
    runner = getattr(_rowmap, 'rowmap', None)

    messages = { 1 : "Computation sucessful.",
                 2 : "Computation interrupted by solout.",
                -1 : "Stepsize too small, i.e. hs < 10*uround*dabs(t).",
                -2 : "More than iwork(1) steps.",
                -3 : "Input is not consistent.",
                -4 : "Internal Krylov matrix is repeatedly singular."
               }

    supports_run_relax = 0
    supports_step = 0
    active_global_handle = 0

    """ Dummy callback functions """
    def dummy_callback(self, *args, **kwargs):
        pass


    """ Simple callback that checks for NaN """
    def check_nan_callback(self, told, tnew, uold, unew, fold, fnew, ucon, intr, ipar):
        if np.any(np.isnan(unew)):
            warnings.warn(self.__class__.__name__ + ': ' + ' Rowmap produced NaN.')
            self.success = 0


    def __init__(self,
                 method='grk4t',
                 with_jacobian=False,
                 rtol=1e-6, atol=1e-12,
                 lbd=0.25, ubd=2., step_safety=0.8,
                 ktol=1.e-1, max_iter=1000, mx=90,
                 lun=6, dt=0.1, *args, **kwargs):

        # set method
        if re.match(method, r'shampine', re.I):
            self.meth = 1
        elif re.match(method, r'grk4t', re.I):
            self.meth = 2
        elif re.match(method, r'gamma', re.I):
            self.meth = 3
        elif re.match(method, r'D-stable', re.I):
            self.meth = 4
        elif re.match(method, r'L-stable', re.I):
            self.meth = 5
        elif re.match(method, r'grk4a', re.I):
            self.meth = 6
        else:
            raise ValueError('Unknown integration method %s' % method)

        self.with_jacobian  = with_jacobian
        self.rtol           = rtol
        self.atol           = atol

        # lower and upper step bounds
        self.lbd            = lbd
        self.ubd            = ubd
        self.step_safety    = step_safety
        self.ktol           = ktol
        self.max_iter       = max_iter
        self.dt             = dt

        # solver output
        self.lun            = lun

        # max krlov subspace size
        self.mx             = mx
        self.success        = 1

        # sets whether rhs is autonomous or non-autonomous
        self.ifcn           = kwargs.pop('ifcn', 0)

        # external functions
        self.solout         = kwargs.pop('solout', None)
        self.fdt            = kwargs.pop('fdt'   , None)
        self.jacv           = kwargs.pop('jacv'  , None)

        # If we provide a time derivative of f then we must be non-auto
        if self.fdt is not None:
            self.ifcn = 1

        # rowmap control pars map
        self.ccontrol       = {'solout' : 'iout', 'fdt' : 'ifdt', 'jacv' : 'ijacv'}

        # todo these are f, jac etc number of parameters
        self.rpar           = 1.
        self.ipar           = 1

        # debug control
        self.debug          = kwargs.pop('debug', False)
        self.verbose        = kwargs.pop('verbose', False)

        # if debug call NaN check
        if self.debug:
            print('Register NaN check!')
            self.solout     = self.check_nan_callback

        # keep track of init status
        self.initialized    = False


    def _get_iwork(self):
        # prepare iwork see rowmap documentation!
        iwork = np.zeros(self.mx + 20, np.int32)
        assert (self.meth>=1) and (self.meth<=6), ''
        iwork[0] = self.max_iter
        iwork[1] = self.meth
        iwork[2] = self.mx
        iwork[3] = self.lun

        return iwork


    def _get_work(self, n):
        # prepare work
        work = np.zeros(10 + n*(self.mx+11)+self.mx*(self.mx+4))
        work[0] = self.lbd
        work[1] = self.ubd
        work[2] = self.step_safety
        work[3] = np.finfo(float).eps
        work[4] = self.ktol

        return work


    def _compute_callback_options(self, kwargs):
        for callback_name in ['solout', 'fdt', 'jacv']:
            callback = getattr(self, callback_name)
            control_name = self.ccontrol[callback_name]

            kwargs[control_name]  = 1 if callback is not None else 0
            kwargs[callback_name] = callback if callback is not None else self.dummy_callback

        return kwargs


    def reset(self, n, has_jac):
        """ Prepare integrator for call: allocate memory, set flags etc.
        n - number of equations.
        has_jac - user supplied jacobian
        """

        #print('Calling reset n=%d' % n)
        # n - number of equations
        # has_jac - whether user has supplied its own Jacobian
        self.iwork = self._get_iwork()
        self.work  = self._get_work(n)

        self.call_args = [self.dt, self.rtol, self.atol,
                          self.work, self.iwork,
                          self.rpar, self.ipar]

        self.call_kwargs  = {'ifcn' : self.ifcn}

        # check if we should enable any callbacks
        self._compute_callback_options(self.call_kwargs)

        self.success = 1
        self.initialized = False


    def run(self, f, jac, y0, t0, t1, f_params, jac_params):
        # this method is called to integrate from t=t0 to t=t1
        # with initial condition y0. f and jac are used-supplied functions
        # that defined the problem. f_params, jac_params are addition arguments
        # to these functions
        if self.initialized:
            self.check_handle()
        else:
            self.initialized = True
            self.acquire_new_handle()

        args    = ((f, t0, y0, t1) + tuple(self.call_args))

        if self.debug:
            print('Args:', args)
            print('kwargs:', self.call_kwargs)

        t, y1, hs, iwork, idid = self.runner(*args, **self.call_kwargs)

        if idid < 0:
            warnings.warn(self.__class__.__name__ + ': ' +
                          self.messages.get(idid, 'Unexpected idid=%s' % idid))
            self.success = 0

        # safe these values for inspection
        self.hs     = hs
        self.iwork  = iwork
        self.idid   = idid

        # IMPROVE THIS!
        if self.verbose or self.debug:
            print(self.statistics())

        if not self.success:
            print('ROWMAP: Failure!')

        if not self.success or self.debug:
            print(self.rowmap_parameters())
            print(self.statistics())

        return y1, t


    def statistics(self):
        return '\tComputed steps %d\n\tRejected steps %d\n\tFunction evals %d'\
                '\n\tMat-Vec products %d\n\ths %.4g\n\tidid %d.'\
                % (self.iwork[4], self.iwork[5], self.iwork[6], self.iwork[7],
                  self.hs, self.idid)


    def rowmap_parameters(self):
        return '\tRowmap parameters: MaxIter %d\n\tMethod %d\n\tmx %d\n\tlun %d.'\
                % (self.iwork[0], self.iwork[1], self.iwork[2], self.iwork[3])


if rowmap.runner:
    IntegratorBase.integrator_classes.append(rowmap)
