#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 15:17:38 2019

@author: usingh

Read pyrpipe configuration
"""

import json
import os
import multiprocessing
import psutil

#read pyrpipe.conf in current dir

class Conf:
    def __init__( self):
        self._dry = False
        self._safe = False
        self._threads=multiprocessing.cpu_count()
        self._params_dir='./params'
        self._memory=None
        
        conf_file_path='pyrpipe.conf'
        if os.path.exists(conf_file_path):
            with open(conf_file_path) as conf_file:
                data = json.load(conf_file)
            self._dry=data['dry']
            self._threads=data['threads']
            self._params_dir=data['params_dir']
            self._memory=data['memory']
            self._safe = data['safe']
            #check valid threads and mem
            if not self._threads or self._threads.replace('.','',1).isdigit():
                self._threads=multiprocessing.cpu_count()
            if not self._memory or self._memory.replace('.','',1).isdigit():
                self._memory=psutil.virtual_memory()[0]/1000000000
            
                
conf=Conf()
_dryrun=conf._dry
_safe=conf._safe
_threads=str(conf._threads)
_mem=str(conf._memory)
_params_dir=conf._params_dir


