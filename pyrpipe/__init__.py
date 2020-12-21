#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 15:17:38 2019

@author: usingh

Read pyrpipe configuration
"""

import sys
import json
import os
import multiprocessing
import psutil
import platform
from multiprocessing import cpu_count
import logging
from datetime import datetime 
import subprocess
import time
import pyrpipe.version
import pyrpipe.arg_parser


###logger
class LogFormatter():
    """
    A formatter for logs
    """
    def __init__(self):
        self.start_time = time.time()    
    def format(self, record):      
        return "{}".format(record.getMessage())       

class PyrpipeLogger():
    """
    Class to manage pyrpipe logs
    
    Attributes
    -----------    
    env_logger: logger to log the current environment
    cmd_logger: logger to log the execution status, stdout, stderr and runtimes for each command run using execute_command()
    """
    def __init__(self,logdir=None):
        self.__name__="pyrpipeLogger"
        #loggers
        timestamp=str(datetime.now()).replace(" ","-").replace(":","_")
        self.logger_basename=timestamp+"_pyrpipe"
        if not logdir:
            logdir="pyrpipe_logs"
        self.logs_dir=os.path.join(os.getcwd(),logdir)
        if not os.path.isdir(self.logs_dir):
            os.mkdir(self.logs_dir)
        self.logger_path=os.path.join(self.logs_dir,self.logger_basename+".log")
        

        self.log_path=os.path.join(self.logs_dir,self.logger_basename+".log")
        self.envlog_path=os.path.join(self.logs_dir,self.logger_basename+"ENV.log")
        formatter=LogFormatter()
        self.env_logger=self.create_logger("env",self.envlog_path,formatter,logging.DEBUG)
        self.cmd_logger=self.create_logger("cmd",self.log_path,formatter,logging.DEBUG)
        self.init_envlog()
        self.init_cmdlog()        
    
    def create_logger(self,name,logfile,formatter,level=logging.DEBUG):
        """Creates a logger
        
        Parameters
        ----------
        
        name: str
            name of logger
        logfile: str
            file name to save logs
        formatter: formatter object
            formatter for log
        
        Returns: logger
            A logger object
        """
        #Get different loggers
        handler = logging.FileHandler(logfile)        
        handler.setFormatter(formatter)
        
        logger = logging.getLogger(name)
        logger.setLevel(level)
        logger.addHandler(handler)
        return logger
    
    def init_cmdlog(self):
        """init the cmdlog
        """
        self.cmd_logger.debug("# Start Log")    
        self.cmd_logger.debug("# pyrpipe version: "+pyrpipe.version.__version__)    

    def init_envlog(self):
        """init the envlog
        """
        self.env_logger.debug("#Start Log")
        self.env_logger.debug("# pyrpipe version: "+pyrpipe.version.__version__)    
        osInfo=platform.platform()
        #get python version
        pyver='Python ' + sys.version.replace('\n', '')
        #get cpu
        cpu=str(cpu_count())+' logical CPU cores'
        #get current conda environment and dump to log
        condaenv_cmd='conda env export'
        try:
            result = subprocess.Popen(condaenv_cmd.split(),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,)
            stdout,stderr = result.communicate()
            condaenv=stdout.decode("utf-8") 
            #remove prefix line       
            condaenv='\n'.join([ x for x in condaenv.split('\n') if x and 'prefix' not in x ])
            #print(condaenv)
        except:
            condaenv='Conda not found'
        
        envDesc={'now':str(datetime.now().strftime("%y-%m-%d %H:%M:%S")),
                 'python':pyver,
                 'os':osInfo,
                 'cpu':cpu,
                 'syspath':str(sys.path),
                 'sysmodules':str(list(sys.modules.keys())),
                 'conda_env':str(condaenv)
                 }
        
        self.env_logger.debug(json.dumps(envDesc))
        self.env_logger.debug("#PROGRAMS")
        #a list of logged programs
        self.logged_programs=[]

#read pyrpipe.conf in current dir
class Conf:
    def __init__( self):
        #if conf file is not present these default values will be used
        self._dry = False
        self._safe = False
        #roughly 80% of available cpu cores
        self._threads=max(int(multiprocessing.cpu_count()*0.8),1)
        self._force=False
        self._params_dir='params'
        #roughly 80% of available memory
        self._memory=psutil.virtual_memory()[0]/1000000000*0.8
        self._logging=True
        self._logs_dir='pyrpipe_logs'
        self._verbose=False
        
        conf_file_path='pyrpipe.conf'
        if os.path.exists(conf_file_path):
            with open(conf_file_path) as conf_file:
                data = json.load(conf_file)
            self._dry=data['dry']
            self._threads=data['threads']
            self._force=data['force']
            self._params_dir=data['params_dir']
            self._logging=data['logging']
            self._logs_dir=data['logs_dir']
            self._verbose=data['verbose']
            self._memory=data['memory']
            self._safe = data['safe']
            #check valid threads and mem
            if not self._threads or not self._threads.replace('.','',1).isdigit():
                self._threads=max(int(multiprocessing.cpu_count()*0.8),1)
            if not self._memory or not self._memory.replace('.','',1).isdigit():
                self._memory=psutil.virtual_memory()[0]/1000000000*0.8
                
        #TODO: overwrite any arguments paased via cmd line
        #self.init_sys_args()
        
    #Not used    
    def init_sys_args(self):
        args, unknownargs = pyrpipe.arg_parser.parser.parse_known_args() 
        if args.threads:
            self._threads=args.threads
        if args.mem:
            self._memory=args.mem
        self._dry=args.dryrun
        self._safe=args.safemode
        if args.paramdir:
            self._params_dir=args.paramdir
        if args.logsdir:
            self._logs_dir=args.logsdir
        self._logging=not args.nologs
        self._verbose=args.verbose
        self._force=args.force
        
            
                
conf=Conf()
_dryrun=conf._dry
_safe=conf._safe
_threads=str(conf._threads)
_mem=str(conf._memory)
_params_dir=conf._params_dir
_logs_dir=conf._logs_dir
_logging=conf._logging
_verbose=conf._verbose
_force=conf._force

#create logger
#pyrpipe_logger=PyrpipeLogger()
