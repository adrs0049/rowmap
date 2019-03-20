#!/usr/bin/python
# -*- coding: utf-8 -*-
# Author: Andreas Buttenschoen

import sys

class Printer(object):
    """
    Print things to stdout by updating line dynamically
    """
    def __init__(self, data):
        sys.stdout.write("\x1b[2K\r" + data.__str__())
        sys.stdout.flush()

