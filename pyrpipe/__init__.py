#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 15:17:38 2019

@author: usingh

Read pyrpipe configuration
"""

import json
import os

#read pyrpipe.conf in current dir

class Conf:
    def __init__( self):
        self._dry = False
        self._threads=4
        self._params_dir='./params'
        conf_file_path='pyrpipe.conf'
        if os.path.exists(conf_file_path):
            with open(conf_file_path) as conf_file:
                data = json.load(conf_file)
            self._dry=data['dry']
            self._threads=data['threads']
            self._params_dir=data['params_dir']
                
conf=Conf()
_dryrun=conf._dry
_threads=str(conf._threads)
_params_dir=conf._params_dir