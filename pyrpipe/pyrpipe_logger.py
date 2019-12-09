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
    

def createLogger(name,logfile,formatter="",level=logging.DEBUG):

    #Get different loggers
    
    handler = logging.FileHandler(logfile)        
    handler.setFormatter(LogFormatter())
    
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


timestamp=str(datetime.now()).split(".")[0].replace(" ","-").replace(":","_")

#create logs directory
logsDir=os.path.join(os.getcwd(),"pyrpipe_logs")
os.mkdir(logsDir)
#log names
cmdLogFname=os.path.join(logsDir,timestamp+"_pyrpipeCMD.log")
stdOutLogFname=os.path.join(logsDir,timestamp+"_pyrpipeOUT.log")
stdErrLogFname=os.path.join(logsDir,timestamp+"_pyrpipeERR.log")

commandLogger=createLogger("cmd",cmdLogFname)
commandLogger.warn("#START LOG")

stdOutLogger=createLogger("stdout",stdOutLogFname)
stdOutLogger.warn("#START LOG")
                   
stdErrLogger=createLogger("stderr",stdErrLogFname)
stdErrLogger.warn("#START LOG")


print ('\033[94m' + "Logs will be written to {} {} {}:".format(cmdLogFname,stdOutLogFname,stdErrLogFname) + '\033[0m')