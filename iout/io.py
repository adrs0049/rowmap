#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Author: Andreas Buttenschoen
#
from __future__ import absolute_import, print_function, division

from lxml import etree
from iout.config import get_config
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


def load_datafile(datafile):
    assert os.path.exists(datafile), 'File %s does not exist!' % datafile
    groups = get_groups(datafile)

    # if more than a single group in hdf5 file we must split it up!
    dfs = {}
    if len(groups) > 1:
        for group in groups:
            print('datafile:', datafile, ' group:', group)
            dfs[group] = pd.read_hdf(datafile, key=group)
    else:
        dfs = pd.read_hdf(datafile)

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


