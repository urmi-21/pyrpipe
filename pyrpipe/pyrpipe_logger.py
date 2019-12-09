#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 17:16:16 2019

@author: usingh
"""

#for logging

from datetime import datetime 
timestamp=str(datetime.now()).split(".")[0].replace(" ","-").replace(":","_")
import logging
import logaugment


class LogFormatter():
    def __init__(self):
        self.start_time = time.time()    
    

logformat = """%(asctime)s : %(message)s\nApprox time: %(time_diff)s Total time: %(relativeCreated)d"""

#format='Time%(asctime)s: %(time_since_last)s: %(message)s'
    
logFile=timestamp+"_pyrpipe.log"
print("Logs will be saved to:"+logFile)
logging.basicConfig(filename=logFile, 
                    format=logformat,filemode='w') 
logger = logging.getLogger()
logger.setLevel(logging.DEBUG) 

def getApproxTimeFromLog(record):
    now = datetime.utcnow()
    try:
        delta = now - getApproxTimeFromLog.now
    except AttributeError:
        delta = 0
    getApproxTimeFromLog.now = now

    return {'time_diff': delta}

logaugment.add(logger, getApproxTimeFromLog)
logger.warn("#START LOG")


"""
stdOutLogFile=timestamp+"_pyrpipe_STDOUT.log"




stdErrLogFile=timestamp+"_pyrpipe_STDERR.log"

formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')

def createCommandsLogger(name,logfile,level=logging.DEBUG):
    
    handler = logging.FileHandler(logfile)        
    handler.setFormatter(formatter)
    
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger
"""
"""
def createLogger(name,logfile,formatter,level=logging.DEBUG):

    #Get different loggers
    
    handler = logging.FileHandler(logfile)        
    handler.setFormatter(formatter)
    
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger

def getLoggerFormatter():
    return logging.Formatter('%(asctime)s %(levelname)s %(message)s')


commandLogger=createLogger("l1",timestamp+"_pyrpipe.log",getLoggerFormatter())
"""