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
import shutil
import platform
from multiprocessing import cpu_count
import logging
from datetime import datetime 
import subprocess
import time
import pyrpipe.version
import pyrpipe.arg_parser
import atexit 
from pyrpipe import pyrpipe_utils as pu



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
    def __init__(self,name,logdir=None):
        self.__name__="pyrpipeLogger"
        #loggers
        #timestamp=str(datetime.now()).replace(" ","-").replace(":","_")
        self.logger_basename=name
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
        self.cmd_logger.debug("# Script file: "+_scriptfile)    
        self.cmd_logger.debug("# MD5 checksum: "+_md5)    
        #log command line
        self.cmd_logger.debug("# Full command: "+_full_command)    
        #log script commands
        self.cmd_logger.debug("# Script options: "+_script_opts)    
        #write mdf for any input files
        for k,v in _optsmd5.items():
            self.cmd_logger.debug("# Input MD5 checksum {}:{}".format(k,v))    
        
        #TODO

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
        #init default values to use
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
        
        #if conf file is present use it
        conf_file_path='pyrpipe.conf'
        if os.path.exists(conf_file_path):
            pu.print_yellow('Reading from pyrpipe.conf')
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
            self.init_threads_mem()
            
        #else use arguments passed
        else:
            args, unknownargs = pyrpipe.arg_parser.parser.parse_known_args()
            if args.versioninfo:
                ver=pyrpipe.version.__version__
                print("pyrpipe version {}".format(ver))
                sys.exit(0)

            infile=args.infile
            """
            if infile is present implies invoked via pyrpipe command and control will go to the __main__ module
            If infile is absent, it was invoked using python script.py command. In this case all pyrpipe options are read
            and removed from the argv parameters so that passed argv  are available to the script.py
            """
            if not infile:
                #invoked using python file.py [opts]
                #update sys.argv to ignore all pyrpipe specific args
                sys.argv=[sys.argv[0]]+unknownargs
            
            #threads
            self._threads=args.threads
            self._memory=args.mem
            self._dry=args.dryrun
            self._safe=args.safemode
            self._params_dir=args.paramdir
            self._logs_dir=args.logsdir
            self._logging=not args.nologs
            self._verbose=args.verbose
            self._force=args.force
            self.init_threads_mem()
            
                
        #TODO: overwrite any arguments paased via cmd line
        #self.init_sys_args()
    
    def init_threads_mem(self):
        #check valid threads and mem
        if not self._threads or not self._threads.replace('.','',1).isdigit():
            self._threads=max(int(multiprocessing.cpu_count()*0.8),1)
        if not self._memory or not self._memory.replace('.','',1).isdigit():
            self._memory=psutil.virtual_memory()[0]/1000000000*0.8
        
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
        
            
#if pyrpipe_diagnostic is invoked
if sys.argv[0].split('/')[-1]=='pyrpipe_diagnostic':
    #will go to __diagnostic__.main
    pass
elif sys.argv[0].split('/')[-1]=='pyrpipe':
    #will go to __main__.main
    pass
else:
    #this will execute in command was like python <script.py> [opts]
    _full_command='python '+' '.join(sys.argv)
    conf=Conf()
    _dryrun=conf._dry
    _safe=conf._safe
    _threads=str(conf._threads)
    _mem=str(conf._memory)
    _params_dir=conf._params_dir
    _logs_dir=conf._logs_dir
    _timestamp=str(datetime.now()).replace(" ","-").replace(":","_")
    _log_name=_timestamp+"_pyrpipe"
    _logging=conf._logging
    _verbose=conf._verbose
    _force=conf._force
    _scriptfile=sys.argv[0]
    _md5=pu.get_mdf(_scriptfile)
    _configuration_path='.pyrpipe'
    _script_opts=','.join(sys.argv)
    
    #create a copy of the script under .pyrpipe folder
    if not pu.check_paths_exist(_configuration_path):
        pu.mkdir(_configuration_path)
    #name: _pyrpipe_hash_filename
    target='_pyrpipe_'+_md5+'_'+pu.get_filename(_scriptfile)
    target=os.path.join(_configuration_path,target)
    #if file not already exist
    if not pu.check_files_exist(target):
        shutil.copyfile(_scriptfile, target)
        pu.print_yellow('Creating script backup: '+target)
    
    #compute has of any arguments if they are files
    _optsmd5={}
    for s in sys.argv[1:]:
        if pu.check_files_exist(s):
            _optsmd5[s]=pu.get_mdf(s)
    

    @atexit.register 
    def goodbye(): 
        if _dryrun:
            pu.print_yellow("This was a dry run. Logs were saved to {}".format(os.path.join(_logs_dir,_log_name+'.log')))
            return
        pu.print_yellow("Logs were saved to {}".format(os.path.join(_logs_dir,_log_name+'.log')))
        #get summary from log
        
        #export shell commands
        
        #run reports/multiqc if specified
        
    
    
    
    
    