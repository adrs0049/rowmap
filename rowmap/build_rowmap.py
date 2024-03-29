#!/usr/bin/python
# -*- coding: utf-8 -*-
# Author: Andreas Buttenschoen
from __future__ import print_function, division

import subprocess, sys, os


def build_rowmap_module():
    dir_cwd = os.getcwd()
    dir_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(dir_path)

    python_exec = sys.executable
    make_process = subprocess.Popen('make clean && make all python_exec=%s' % python_exec,
                                    stdout=subprocess.PIPE,
                                    shell=True, stderr=subprocess.STDOUT)

    if make_process.wait() != 0:
        err_str = make_process.communicate()[0]
        raise NameError('Build of rowmap module failed!\n{}'.format(err_str))

    # switch back to old dir
    os.chdir(dir_cwd)


