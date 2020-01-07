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
from iout.format_time import format_delta, now

import os, time
import numpy as np
import pandas as pd
from scipy.integrate import ode
import rowmap.rowmap_ode_runner

from iout.io import writeDataFrame
from iout.h5_io import MOLFile
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
        self.miter  = kwargs.pop('max_iter', 1000)
        self.tout   = np.arange(self.time.t0, self.time.tf, self.time.dt)

        # lambdas to create initial condition
        self.y0     = y0

        # make sure y0 has the correct shape
        self._reshape_y0()

        # look up domain and dimensions
        self.dim    = self._get_dims(*args, **kwargs)

        # no-equations
        self.noEqs = self._get_no_eqns()

        # live plotting
        self.livePlotting = kwargs.pop('livePlotting', False)
        self.plotter      = None

        # data output
        self.outdir         = kwargs.pop('outdir', 'results')
        self.save           = kwargs.pop('save', False)
        self.save_new       = kwargs.pop('save_new', False)
        assert (self.save and self.save_new) is False, 'Save and new_save cannot be used at the same time!'
        self.name           = kwargs.pop('name'  , 'MOL_unnamed')

        # set default verbosity
        self.verbose        = kwargs.pop('verbose', False)
        self.debug          = kwargs.pop('debug',   False)

        # integrator
        self.ode    = None

        # The right hand side
        self.f      = f(noPDEs = self.noEqs, t0 = self.time.t0, *args, **kwargs)

        # number of spatial disc points
        self.msize = self._get_spatial_pts()

        # if not h5 output save a dataframe
        self._setup_output(y0)

        # machine eps
        self.eps    = 1.e4 * np.finfo(float).eps

        # current solution
        self.yt = None

        # output counter
        self.output_counter = 0
        self.output_every   = kwargs.pop('output_every', 20)

        # setup
        self._setup()


    """ Easy access to current solution """
    @property
    def y(self):
        return self.f.reshape(self.yt)


    """ reshape y0 if it's a 1D array """
    def _reshape_y0(self):
        shape = self.y0.shape
        if len(shape) == 1:
            self.y0 = self.y0.reshape((1, self.y0.size))


    """ Determine the problems dimension """
    def _get_dims(self, *args, **kwargs):
        domain = kwargs.get('domain', None)
        if domain is None:
            # TODO proper warnings!
            print('WARNING: couldn\'t find a domain!')
            return 1

        return domain.dimensions


    """ Determine the number of equations """
    def _get_no_eqns(self):
        shape = self.y0.shape

        if self.dim > 2:
            assert False, 'Dimensions other than 1 and 2 are not support!'
        else:
            if len(shape) > self.dim:
                return shape[0]
            else:
                # otherwise we must have 1 PDE only!
                return 1


    """ Determine number of spatial points """
    def _get_spatial_pts(self):
        shape = self.y0.shape

        if self.dim > 2:
            assert False, 'Dimensions other than 1 and 2 are not support!'
        else:
            return shape[1:]


    """ setup dataframe output """
    def _create_outdir(self):
        if os.path.exists(self.outdir) and not os.path.isdir(self.outdir):
            assert False, '%s exists but is not a directory!' % self.outdir

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)


    """ setup output """
    def _setup_output(self, y0):
        if self.save_new:
            self.write_disk = self.write_disk_h5
            self.write = self.write_h5
            self._setup_h5_output(y0)
        else:
            self.write_disk = self.write_disk_df
            self.write= self._writeToLocalDataFrame
            self._setup_df_output(y0)


    """ Use a pandas dataframe for outputs """
    def _setup_df_output(self, y0):
        self._create_outdir()

        # data storage
        self.dfs = {}
        if self.noEqs == 1:
            y0 = np.reshape(y0, (1, y0.size))

        for i in range(self.noEqs):
            y = y0[i, :].flatten()
            # creating it with columns is super slow
            self.dfs[i] = pd.DataFrame(index=range(y.size)).transpose()
            self.dfs[i].name = 'MOL_dataframe' + str(i)
            self.dfs[i].loc[float(self.time.t0)] = y

        # write to HDF5-file
        self.write_disk(append=False)


    """ Direct HDF5 output """
    def _setup_h5_output(self, y0):
        self._create_outdir()

        self.h5f = MOLFile(shape=self.msize, eqns=self.noEqs, fname=self.outfile,
                           name=self.name)

        if self.noEqs == 1:
            y0 = np.reshape(y0, (1, y0.size))

        # save initial condition to hdf5 file
        self.h5f.add(float(self.time.t0), y0)


    """ update the plotter """
    def update_plot(self, t, y):
        # update the plot
        if self.livePlotting:
            self.plotter.update(self.f.dom.xs(), y, title = 't = %.2f' % t)


    """ print info """
    def __str__(self):
        return 'MOL(%s).' % (self.time)


    """ Get the result HDF5 file """
    @property
    def outfile(self):
        return os.path.join(self.outdir, self.name + '.h5')


    """ Write local dataframe to a HDF5 file """
    def write_disk_df(self, append=True):
        if not self.save:
            return

        # TODO make this general!
        for df in self.dfs.values():
            writeDataFrame(self.outfile, df, append=append)

        # now drop all data from dataframe
        for df in self.dfs.values():
            df.drop(df.index, inplace=True)


    """ Write to h5 file """
    def write_h5(self, y, t):
        self.h5f.add(t, y)


    """ TODO remove! """
    def write_disk_h5(self, append=True):
        pass


    """ Update the local dataframe """
    def _writeToLocalDataFrame(self, y, t):
        # write all data to the local dataframe
        for i in range(self.noEqs):
            yy = y[i, :].flatten()
            self.dfs[i].loc[t] = yy


    """ Integrate using MOL """
    def run(self):
        start       = now()
        step_start  = now()

        #print("MOL TIME:", self.time)
        while self.ode.successful() and self.time.keepGoing(self.ode.t):
            try:
                self.yt = self.ode.integrate(self.ode.t + self.time.dt)

            # Handle value errors differently
            except ValueError as e:
                # force a write of the current state
                if self.yt is not None:
                    self.write_disk()
                    yf = self.f.reshape(self.yt)
                    self.write(yf, self.ode.t)

                    # Write to HDF5-file
                    self.write_disk()
                else:
                    print('MOL error: Solver returned None type!\n\tCan\'t save state!')

                raise e

            # catch anything else
            except Exception as e:
                raise e

            # reformat yt
            yf = self.f.reshape(self.yt)

            self.write(yf, self.ode.t)

            # Check if we should write to the disk
            self.output_counter += 1

            if (self.output_counter % self.output_every) == 0:
                self.output_counter = 0

                # Write to HDF5-file
                self.write_disk()

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
        self.write_disk()

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
        # check if rhs is autonomous
        ifcn = 0 if self.f.isAutonomous() else 1

        self.ode = ode(self.f).set_integrator('rowmap', method='grk4t', dt=self.hi,
                                              rtol=self.vtol, ifcn=ifcn,
                                              atol=self.vtol**2,
                                              ktol = self.ktol,
                                              max_iter=self.miter,
                                              debug=self.debug)

        self.ode.set_initial_value(self.y0.flatten(), self.time.t0)


if __name__ == '__main__':
    rowmap = rowmap.rowmap_ode_runner
    print(rowmap.__version__)

