#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Author: Andreas Buttenschoen
#
from __future__ import absolute_import, print_function, division

from lxml import etree
from iout.config import get_config
from iout.h5_io import MOLFile

from collections import OrderedDict

import os
import h5py as h5
import pandas as pd


def get_groups(datafile):
    assert os.path.exists(datafile), 'File %s does not exist!' % datafile
    f = h5.File(datafile, 'r')
    groups = list(f.keys())
    close_h5(f)
    return groups


def close_h5(fhandler):
    if isinstance(fhandler, h5.File):
        try:
            fhandler.close()
        except:
            pass # was already closed


""" Function to load a hdf5 into pytables """
def load_datafile_pytables(datafile):
    groups = get_groups(datafile)

    # if more than a single group in hdf5 file we must split it up!
    dfs = {}
    if len(groups) > 1:
        for group in groups:
            dfs[group] = pd.read_hdf(datafile, key=group)
    else:
        dfs = pd.read_hdf(datafile)

    return dfs


""" Load custom hdf5 files """
def load_datafile_mol(datafile):
    print('Loading MOL file: %s.' % datafile)
    mfn = MOLFile(fname=datafile, fmode='r')

    # get keys
    keys = set(mfn.keys())

    # check that we found enough keys
    assert len(keys) > 1, 'Only found %d keys, not enough! Need at least 2.' % len(keys)
    assert 'time' in keys, 'Could not find a group containing simulation time!'

    # load simulation times
    times = mfn['time']

    # if we only have one field we don't need to generate a dictionary!
    fields = keys.symmetric_difference(['time'])

    # the output dictionary
    dfs = OrderedDict()

    if len(fields)>1:
        for group in fields:
            data       = mfn[group+'/data']
            dfs[group] = pd.DataFrame(index=times, data=data)
    else:
        data = mfn[next(iter(fields))+'/data']
        dfs = pd.DataFrame(index=times, data=data)

    return dfs


def load_datafile(datafile):
    assert os.path.exists(datafile), 'File %s does not exist!' % datafile

    # first try to use pytables
    try:
        dfs = load_datafile_pytables(datafile)
    except Exception as e:
        # try to load using molfile
        try:
            dfs = load_datafile_mol(datafile)
        except Exception as e:
            print('Loading %s failed both using MOL and pytables!' % datafile)
            raise e

    return dfs


def writeDataFrame(filename, dataframe, *args, **kwargs):
    # writing a dataframe to hdf5
    dataframe.to_hdf(filename, dataframe.name, mode='a', table=True,
                     format='table', *args, **kwargs)


def print_attrs(name, obj):
    print(name)
    for key, val in obj.attrs.iteritems():
        print("    %s: %s" % (key, val))


def getXMLNode(tree, name):
    root = tree.getroot()
    for child in root:
        if child.tag == name:
            return child


def parseXML(xmlFile, name='Simulator'):
    assert os.path.isfile(xmlFile), '%s not found!' % xmlFile
    parser = etree.XMLParser(remove_comments=True)
    tree   = etree.parse(xmlFile, parser)
    return getXMLNode(tree, name)


def setup_project(projectName):
    config = get_config(projectName)
    create_dir(config['DEFAULT']['outdir'])


def create_dir(dirname):
    dirname=os.path.expanduser(dirname)
    if not os.path.exists(dirname):
        os.makedirs(dirname)


def openH5File(dataFile, verbose=False):
    if verbose:
        print('Trying to open %s' % dataFile)

    assert os.path.isfile(dataFile), 'Datafile %s not found.' % dataFile
    fn = h5.File(dataFile, 'r')
    outputPath = os.path.dirname(os.path.abspath(dataFile))

    if verbose:
        fn.visititems(print_attrs)
        print('Writing results to %s.' % outputPath)

    return fn, outputPath
