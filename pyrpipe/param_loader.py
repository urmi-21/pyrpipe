#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 14:37:07 2020

@author: usingh
"""

import yaml
from pyrpipe import pyrpipe_utils as pu

#ignore boolean coversion of yes/no values when loading yaml
from yaml.constructor import SafeConstructor
def add_bool(self, node):
    return self.construct_scalar(node)
SafeConstructor.add_constructor(u'tag:yaml.org,2002:bool', add_bool)

class YAML_loader():
    """
    Load parameters from a yaml file
    """
    
    
    def __init__(self,file):
        self.__params=None
        self.__kwargs=None
        
        if not pu.check_files_exist(file):
            return
        #read yaml
        with open(file) as f:
            #self.__params=yaml.full_load(f)
            self.__params= yaml.load(f, Loader=yaml.SafeLoader)
            self.parse_params()
        

    def get_params(self):
        return self.__params
    
    def get_kwargs(self):
        if self.__kwargs:
            return self.__kwargs
        return {}
    
    def parse_params(self):
        """
        store params as dict
        """
        #if file is empty
        if not self.__params:
            self.__kwargs={}
            return
        #create copy
        params=self.__params.copy()
        to_del=[]
        for k,v in params.items():
            #handle boolean arguments
            if type(v)==type(True) and v == True:
                params[k]=""
            elif type(v)==type(True) and v == False:
                to_del.append(k)
            else:
                #convert int, num to string
                params[k]=str(v)
        
        for k in to_del:
            del params[k]        
        
        self.__kwargs=params
        
        
                