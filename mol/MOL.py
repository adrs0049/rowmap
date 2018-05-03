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

import os, time
import numpy as np
import pandas as pd
from scipy.integrate import ode
import rowmap.rowmap_ode_runner

from iout.io import writeDataFrame
from iout.Printer import Printer

from model.Time import Time

# visualization
from vis.plot import Plot


# TODO link these parameter classes directly to the mysql models

class MOL:
    def __init__(self, f, y0, *args, **kwargs):
        self.version    = 'MOL-0.1'

        # time control
        self.time   = kwargs.pop('time', Time())

        # integrator control
        self.hi     = kwargs.pop('hi', 1e-3)
        self.vtol   = kwargs.pop('vtol', 1e-3)
        self.ktol   = kwargs.pop('ktol', 1e-1)
        self.tout   = np.arange(self.time.t0, self.time.tf, self.time.dt)

        # lambdas to create initial condition
        self.y0     = y0

        # no-equations
        self.noEqs = self._get_no_eqns()

        # data storage
        self.dfs = {}
        if self.noEqs == 1:
            y0 = np.reshape(y0, (1, y0.size))

        for i in range(self.noEqs):
            self.dfs[i] = pd.DataFrame()
            self.dfs[i].name = 'MOL_dataframe' + str(i)
            self.dfs[i][float(self.time.t0)] = y0[i, :].flatten()

        # live plotting
        self.livePlotting = kwargs.pop('livePlotting', False)
        self.plotter      = None

        # data output
        self.outdir         = kwargs.pop('outdir', 'results')
        self.name           = kwargs.pop('name'  , 'MOL_unnamed')
        self.save           = kwargs.pop('save', False)

        # set default verbosity
        self.verbose        = kwargs.pop('verbose', False)

        # integrator
        self.ode    = None


        # The right hand side
        self.f      = f(*args, **kwargs)

        # machine eps
        self.eps    = 1.e4 * np.finfo(float).eps

        # current solution
        self.yt = None

        # setup
        self._setup()


    """ Determine the number of equations """
    def _get_no_eqns(self):
        shape = self.y0.shape
        print('shape:',shape)
        if len(shape) > 1:
            return shape[0]
        else:
            return 1


    """ update the plotter """
    def update_plot(self, t, y):
        # update the plot
        if self.livePlotting:
            self.plotter.update(self.f.dom.xs(), y, title = 't = %.2f' % t)


    """ print info """
    def __str__(self):
        return 'MOL(%s).' % (self.time)


    def write(self):
        if self.save:
            writeDataFrame(os.path.join(self.outdir, self.name + '.h5'), self.df)


    """ Update the local dataframe """
    def _writeToLocalDataFrame(self, y, t):
        for i in range(self.noEqs):
            yy = y[i, :].flatten()
            self.dfs[i][t] = yy


    """ Integrate using MOL """
    def run(self):
        start       = now()
        step_start  = now()

        #print("MOL TIME:", self.time)
        while self.ode.successful() and self.time.keepGoing(self.ode.t):
            self.yt = self.ode.integrate(self.ode.t + self.time.dt)

            # reformat yt
            yf = self.f.reshape(self.yt)

            self._writeToLocalDataFrame(yf, self.ode.t)

            # tell the runner something about the status
            #if self.ode.t > iteration * tenth_of_run_time:
            end = now()
            if self.verbose:
                ostr = 'Simulation time: %.2f of %.2f in %s (step %s).' \
                    % (self.ode.t, self.time.tf, format_delta(start, end),
                       format_delta(step_start, end))

            # update the plot
            self.update_plot(self.ode.t, yf)

            if self.verbose: Printer(ostr)
            step_start = now()

        # write everything now
        self.write()

        # print new line
        #print('\n')


    def solution(self):
        return self.f.reshape(self.yt)


    #""" Change sizes """
    #def resize(self, newsize):
    #    # newsize tells use how many
    #    if newsize > self.y0.shape[1]:
    #        # Have to grow stuff
    #        self.f.resizeDomain()
    #        self.y0 = np.lib.pad(self.y0, 'constant', constant_values=(0))
    #
    #    elif newsize < self.y0.shape[1]:
    #        # have to shrink stuff
    #        self.f.resizeDomain()
    #
    #    self._setup()
    #    # UPDATE self.f


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

        if self.livePlotting:
            self.plotter = Plot(self.f.dom.box(), labels=['a', 'b'])
            self.plotter.initialize(self.f.dom.xs(), self.f.reshape(self.y0), title = "t = 0")
            time.sleep(1.)


    def _setup_integrator(self):
        self.ode = ode(self.f).set_integrator('rowmap', method='grk4t', dt=self.hi,
                                              rtol=self.vtol,
                                              atol=self.vtol**2,
                                              ktol = self.ktol)

        self.ode.set_initial_value(self.y0.flatten(), self.time.t0)

