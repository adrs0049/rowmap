#!/usr/bin/python
# -*- coding: utf-8 -*-
# Author: Andreas Buttenschoen

import numpy as np


class Time(object):
    def __init__(self, *args, **kwargs):
        self.tf = kwargs.pop('tf', 1.)
        self.dt = kwargs.pop('dt', 0.1)
        self.t0 = kwargs.pop('t0', 0.)

        self.currentTime = self.t0

        # machine eps
        self.eps    = 1.e4 * np.finfo(float).eps


    def reset(self):
        self.currentTime = self.t0


    def step(self):
        self.currentTime += self.dt


    def keepGoing(self, t):
        return np.abs(self.tf - t) > self.eps


    def __call__(self):
        return self.currentTime


    def __str__(self):
        return 'Time(tf = %.2f, dt = %.2f, t0 = %.2f)' % (self.tf, self.dt, self.t0)


    def __repr__(self):
        return self.__str__()


    def __iadd__(self, other):
        self.tf += other
        self.t0 += other
        return self

