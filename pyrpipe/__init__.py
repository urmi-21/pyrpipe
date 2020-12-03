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
        conf_file_path='pyrpipe.conf'
        if os.path.exists(conf_file_path):
            with open(conf_file_path) as conf_file:
                data = json.load(conf_file)
            self._dry=data['dry']
        
    #read conf file
    @classmethod
    def set( cls, value, label = 'default' ):
        with open('pyrpipe.conf') as conf_file:
            data = json.load(conf_file)
        cls._dry=data['dry']
        
dryrun=Conf()._dry