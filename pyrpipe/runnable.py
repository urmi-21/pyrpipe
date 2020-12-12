#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 13:55:13 2020

@author: usingh
"""
import os
from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe
from pyrpipe import param_loader as pl
from pyrpipe import _params_dir

class Runnable:
    """The runnable class
    """
    def __init__(self,*args,_command=None,**kwargs):
        self._runnable=True
        self._command=None
        self._deps=None
        self._param_yaml=None
        self.init_parameters(*args,**kwargs)
        self._args_style='LINUX'
        self._valid_args=None
         
    def check_dependency(self):
        if self._deps:
            if not pe.check_dependencies(self._deps):
                raise Exception("ERROR. Please check dependencies for {}. Deps: {}".format(self._command," ".join(self._deps)))    
        
    def init_parameters(self,*args,**kwargs):
        #init the parameters for the object
        if args:
            self._args=args
        else:
            self._args=()
        if kwargs:
            self._kwargs=kwargs
        else:
            self._kwargs={}
        #read yaml parameters
        if not self._param_yaml:
            return
        yamlfile=os.path.join(_params_dir,self._param_yaml)
        #override yaml parameters by kwargs
        if pu.check_files_exist(yamlfile):
            yaml_params=pl.YAML_loader(yamlfile)
            yaml_kwargs=yaml_params.get_kwargs()
            self._kwargs={**yaml_kwargs,**self._kwargs}
        
    def run(self,*args,objectid="NA",**kwargs):
        #override class kwargs by passed kwargs
        kwargs={**self._kwargs,**kwargs}
        #if no args provided use constructor's args
        if not args:
            args=self._args
        #if args are not None
        if args and args[0]:
            kwargs['--']=args
        
        if isinstance(self._command, str):
            cmd=[self._command]
        else:
            cmd=self._command
            
        #add options
        if self._args_style=='LINUX':
            cmd.extend(pu.parse_unix_args(self._valid_args,kwargs))      
        elif self._args_style=='JAVA':
            cmd.extend(pu.parse_java_args(self._valid_args,kwargs))
        else: 
            pu.print_boldred("Unknown args style: {}".format(self._args_style))
            raise Exception("Unknown args style")
            
            
        #execute command
        cmd_status=pe.execute_command(cmd,objectid=objectid)
        if not cmd_status:
            pu.print_boldred("{} failed: {}".format(self.name," ".join(cmd)))
     
        #return status
        return cmd_status