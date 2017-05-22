#!/usr/bin/python
# -*- coding: utf-8 -*-
# Author: Andreas Buttenschoen
"""
    This class implements solving PDEs using the Method of Lines.

    The MOL ODE is given by

     d y(t)[i]
    ---------  = f(t, y(t))[i].
       d t

	y(t=0)[i] = y0[i],

	where::

    i = 0, ..., len(y0) - 1

    This ODE is currently solved using the ROWMAP algorithm. For details see
    its implementation in rowmap_ode_runner.py.

    The main object this class expects is a class that implements f(t, y(t)).
    This object is required to look like this.

    class f:
        def __init__(self, *args, **kwargs):
            pass

        def __call__(self, t, y, *args, **kwargs):
            # compute new ydot
            return ydot

    Then the main program looks like.

    f = f(...)
    solver = MOL(f)
    solver.run()

"""
from __future__ import print_function, division

import numpy as np
import pandas as pd
from scipy.integrate import ode
import rowmap.rowmap_ode_runner

# TODO link these parameter classes directly to the mysql models

class MOL:
    def __init__(self, f, y0, *args, **kwargs):
        self.version    = 'MOL-0.1'

        # integrator control
        self.t0     = kwargs.pop('t0', 0.)
        self.tf     = kwargs.pop('tf', 1.)
        self.dt     = kwargs.pop('dt', 0.1)
        self.tout   = np.arange(self.t0, self.tf, self.dt)

        # lambdas to create initial condition
        self.y0     = y0

        # data storage
        self.df     = pd.DataFrame()

        # integrator
        self.ode    = None

        # The right hand side
        self.f      = f(*args, **kwargs)

        # setup
        self._setup()


    """ Integrate using MOL """
    def run(self):
        while self.ode.successful() and self.ode.t < self.tf:
            self.df[self.ode.t] = self.ode.integrate(self.ode.t + self.dt)


    """ Internals """
    def _create_mol(self):
        for patchId in self.grd.patches.keys():
            # todo per patch
            self.data.set_values(self.y)
            self.tdr.setup(patchId)


    def _setup(self):
        self._setup_integrator()
        # set ic in tdr
        self.f.update(self.t0, self.y0)


    def _setup_integrator(self):
        self.ode = ode(self.f).set_integrator('rowmap', method='grk4t')
        self.ode.set_initial_value(self.y0, self.t0)




