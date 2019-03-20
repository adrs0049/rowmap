#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Author: Andreas Buttenschoen
#
from __future__ import print_function, division

import configparser
import os


def get_config(projectname, prefix=None):
    config_file = None

    if prefix is None:
        homedir = os.path.expanduser('~')
        config_file = os.path.join(homedir, 'sources', 'NumericalPDEs', 'config', projectname+'.ini')
        assert os.path.isfile(config_file), '%s does not exist!' % config_file
    else:
        assert False, 'Not implemented!'

    return open_config(config_file)


def open_config(filename):
    config = configparser.ConfigParser(allow_no_value=True)
    config.read(filename)
    return config


def print_config(config):
    for section in config.sections():
        print('Section: %s' % section)
        for options in config.options(section):
            print('x %s::%s::%s' % (options, config.get(section, options), str(type(options))))
