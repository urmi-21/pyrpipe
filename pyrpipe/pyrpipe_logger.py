#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 17:16:16 2019

@author: usingh
"""

#for logging
import time
from datetime import datetime 
from datetime import timedelta
import logging
import os
from sinfo import sinfo
import platform
from multiprocessing import cpu_count

##for sinfo output
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import sys

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout

class LogFormatter():
    def __init__(self):
        self.start_time = time.time()
    
    def format(self, record):
        timeNow=str(datetime.now())
        elapsed_seconds = record.created - self.start_time
        hours, remainder = divmod(elapsed_seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        elapsedTime = timedelta(seconds = elapsed_seconds)
        
        return "Time:{} \n{} \nDuration: {:02}:{:02}:{:02}".format(timeNow, record.getMessage(), int(hours), int(minutes), int(seconds) )    
    

def createLogger(name,logfile,formatter,level=logging.DEBUG):
    #Get different loggers
    handler = logging.FileHandler(logfile)        
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


timestamp=str(datetime.now()).split(".")[0].replace(" ","-").replace(":","_")

#create logs directory

logsDir=os.path.join(os.getcwd(),"pyrpipe_logs")
if not os.path.isdir(logsDir):
    os.mkdir(logsDir)
#log names
cmdLogFname=os.path.join(logsDir,timestamp+"_pyrpipeCMD.log")
stdOutLogFname=os.path.join(logsDir,timestamp+"_pyrpipeOUT.log")
stdErrLogFname=os.path.join(logsDir,timestamp+"_pyrpipeERR.log")
programLogFname=os.path.join(logsDir,timestamp+"_pyrpipeEnv.log")

commandLogger=createLogger("cmd",cmdLogFname,LogFormatter())
commandLogger.debug("#START LOG")

stdOutLogger=createLogger("stdout",stdOutLogFname,LogFormatter())
stdOutLogger.debug("#START LOG")
                   
stdErrLogger=createLogger("stderr",stdErrLogFname,LogFormatter())
stdErrLogger.debug("#START LOG")

with Capturing() as output:
    sinfo(write_req_file=False)
envLogger=createLogger("env",programLogFname,logging.Formatter("%(message)s"))
envLogger.debug("#START LOG")

envLogger.debug("\n".join(output))
envLogger.debug("sys.path:"+str(sys.path))
envLogger.debug("#PROGRAMS")
#a list of logged programs
loggedPrograms=[]



print ('\033[94m' + "Logs will be written to {}, {}, {}, {}".format(cmdLogFname,stdOutLogFname,stdErrLogFname,programLogFname) + '\033[0m')