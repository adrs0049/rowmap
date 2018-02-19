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
from .format_time import format_delta, now

import os
import numpy as np
import pandas as pd
from scipy.integrate import ode
import rowmap.rowmap_ode_runner
from iout.io import writeDataFrame
from iout.Printer import Printer

# TODO link these parameter classes directly to the mysql models

class MOL:
    def __init__(self, f, y0, *args, **kwargs):
        self.version    = 'MOL-0.1'

        # integrator control
        self.t0     = kwargs.pop('t0', 0.)
        self.tf     = kwargs.pop('tf', 1.)
        self.dt     = kwargs.pop('dt', 1e-3)
        self.hi     = kwargs.pop('hi', 1e-3)
        self.vtol   = kwargs.pop('vtol', 1e-3)
        self.tout   = np.arange(self.t0, self.tf, self.dt)

        # lambdas to create initial condition
        self.y0     = y0

        # data storage
        self.df             = pd.DataFrame()
        self.df.name        = 'MOL_dataframe'
        # save initial condition as well
        self.df[self.t0]    = y0

        # data output
        self.outdir         = kwargs.pop('outdir', 'results')
        self.name           = kwargs.pop('name'  , 'MOL_unnamed')

        # integrator
        self.ode    = None

        # The right hand side
        self.f      = f(*args, **kwargs)

        # setup
        self._setup()


    """ print info """
    def __str__(self):
        return 'MOL(t0 = %.2g, t = %.2g, tf = %.2g).' % (self.t0, self.ode.t, self.tf)


    def write(self):
        writeDataFrame(os.path.join(self.outdir, self.name + '.h5'), self.df)


    """ Integrate using MOL """
    def run(self):
        start       = now()
        step_start  = now()

        while self.ode.successful() and self.ode.t <= self.tf:
            self.df[self.ode.t] = self.ode.integrate(self.ode.t + self.dt)

            # tell the runner something about the status
            #if self.ode.t > iteration * tenth_of_run_time:
            end = now()
            ostr = 'Simulation time: %.2f of %.2f in %s (step %s).' \
                  % (self.ode.t, self.tf, format_delta(start, end),
                     format_delta(step_start, end))

            Printer(ostr)
            step_start = now()

        # write everything now
        self.write()


    """ Internals """
    def _create_mol(self):
        for patchId in self.grd.patches.keys():
            # todo per patch
            self.data.set_values(self.y)
            self.tdr.setup(patchId)


    def _setup(self):
        self._setup_integrator()
        # set ic in tdr
        #self.f.update(self.t0, self.y0)


    def _setup_integrator(self):
        self.ode = ode(self.f).set_integrator('rowmap', method='grk4t',
                                              dt=self.hi,
                                              rtol=self.vtol, atol=self.vtol**2)
        self.ode.set_initial_value(self.y0.flatten(), self.t0)




