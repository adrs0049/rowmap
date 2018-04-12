#!/usr/bin/python
# -*- coding: utf-8 -*-

#import h5py
#import numpy as np


def writeDataFrame(filename, dataframe):
    # writing a dataframe to hdf5
    dataframe.to_hdf(filename, dataframe.name, mode='a', table=True)


def print_attrs(name, obj):
    print(name)
    for key, val in obj.attrs.iteritems():
        print("    %s: %s" % (key, val))
