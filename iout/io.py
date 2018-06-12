#!/usr/bin/python
# -*- coding: utf-8 -*-

from lxml import etree
#import h5py
#import numpy as np

def writeDataFrame(filename, dataframe):
    # writing a dataframe to hdf5
    dataframe.to_hdf(filename, dataframe.name, mode='a', table=True)


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
    parser = etree.XMLParser(remove_comments=True)
    tree   = etree.parse(xmlFile, parser)
    return getXMLNode(tree, name)

