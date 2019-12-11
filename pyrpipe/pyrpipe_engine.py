#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 16:26:02 2019

@author: usingh

Methods and classes related to execution and logging
"""

import os
import subprocess
import time
from datetime import datetime 
import logging
import sys
import platform
from multiprocessing import cpu_count
from pyrpipe.pyrpipe_utils import *

class LogFormatter():
    def __init__(self):
        self.start_time = time.time()
    
    def format(self, record):
        
               
        timeNow=str(datetime.now())
        elapsed_seconds = record.created - self.start_time
        hours, remainder = divmod(elapsed_seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        
        message=record.getMessage()
        if message.startswith("#"):
            return "{}\t{}".format(record.getMessage(),timeNow)
        
        return "Time:{} \n{} \nDuration: {:02}:{:02}:{:02}".format(timeNow, record.getMessage(), int(hours), int(minutes), int(seconds) )

class pyrpipeLogger():
    def __init__(self):
        self.__name__="pyrpipeLogger"
        #loggers
        timestamp=str(datetime.now()).split(".")[0].replace(" ","-").replace(":","_")
        self.loggerBaseName=timestamp+"_pyrpipe."
        self.logsDir=os.path.join(os.getcwd(),"pyrpipe_logs")
        if not os.path.isdir(self.logsDir):
            os.mkdir(self.logsDir)
            
        self.cmdLoggerPath=os.path.join(self.logsDir,self.loggerBaseName+"CMD.log")
        self.envLoggerPath=os.path.join(self.logsDir,self.loggerBaseName+"ENV.log")
        self.stdoutLoggerPath=os.path.join(self.logsDir,self.loggerBaseName+"OUT.log")
        self.stderrLoggerPath=os.path.join(self.logsDir,self.loggerBaseName+"ERR.log")
        
        
        
        self.cmdLogger=self.createLogger("cmd",self.cmdLoggerPath,LogFormatter(),logging.DEBUG)
        self.envLogger=self.createLogger("env",self.envLoggerPath,logging.Formatter("%(message)s"),logging.DEBUG)
        self.stdoutLogger=self.createLogger("out",self.stdoutLoggerPath,LogFormatter(),logging.DEBUG)
        self.stderrLogger=self.createLogger("err",self.stderrLoggerPath,LogFormatter(),logging.DEBUG)
        
        self.initCmdlog()
        self.initOutlog()
        self.initErrlog()
        self.initEnvlog()
    
    def createLogger(self,name,logfile,formatter,level=logging.DEBUG):
        #Get different loggers
        handler = logging.FileHandler(logfile)        
        handler.setFormatter(formatter)
        
        logger = logging.getLogger(name)
        logger.setLevel(level)
        logger.addHandler(handler)
        return logger
    
    def initCmdlog(self):
        self.cmdLogger.debug("#START LOG")
    def initOutlog(self):
        self.stdoutLogger.debug("#START LOG")
    def initErrlog(self):
        self.stderrLogger.debug("#START LOG")
    def initEnvlog(self):
        self.envLogger.debug("#START LOG")      
        #get current time
        sesstime='Session information collected on {}'.format(datetime.now().strftime('%Y-%m-%d %H:%M'))
        #get os
        osInfo=platform.platform()
        #get python version
        pyver='Python ' + sys.version.replace('\n', '')
        #get cpu
        cpu=str(cpu_count())+' logical CPU cores'

        self.envLogger.debug(sesstime)
        self.envLogger.debug(pyver)
        self.envLogger.debug(osInfo)
        self.envLogger.debug(cpu)
        self.envLogger.debug("#SYS PATH")
        self.envLogger.debug("sys.path:"+str(sys.path))
        self.envLogger.debug("#SYS MODULES")
        self.envLogger.debug("sys.modules:"+str(sys.modules.keys()))
        self.envLogger.debug("#PROGRAMS")
        #a list of logged programs
        self.loggedPrograms=[]
        

###create logger
pyrpipeLoggerObject=pyrpipeLogger()
printYellow("Logs will be saved to {}*.log".format(pyrpipeLoggerObject.loggerBaseName))
    
"""
All functions that interact with shell are defined here. 
"""
    
##############Functions###########################
    


def getCommandReturnStatus(cmd):
    #not logging these commands
    return executeCommand(cmd,logs=False)

#prints stdout in real time. optimal for huge stdout and no stderr
def executeCommandRealtime(cmd):
    pl.logger.debug("Executing command:\n$ q"+" ".join(cmd)) 
    start_time = time.time()
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    
    for stdout_line in iter(popen.stdout.readline, ""):
            yield stdout_line 
    popen.stdout.close()
    return_code = popen.wait()
    end_time = time.time()
    pl.logger.debug("Executing command:\n$ q"+str(end_time - start_time)) 
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


def executeCommand(cmd,verbose=False,quiet=False,logs=True):
    """
    Function to execute commands using popen. All logs are managed inside the function for all the commands executed.
    
    Parameters
    ----------
    cmd: list
        command to execute in a list
    verbose
        whether to print stdout and stderr. Default: False. All stdout and stderr will be saved to logs regardless of this flag.
    """
    logMessage="$ "+" ".join(cmd)
    if not quiet:
        printBlue(logMessage)
    timeStart = time.time()
    try:
        result = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        stdout,stderr = result.communicate()
        #convert to string
        if stdout:
            stdout=stdout.decode("utf-8")
        else:
            stdout=""
        if stderr:
            stderr=stderr.decode("utf-8")
        else:
            stderr=""
        
        timeDiff = time.time() - timeStart
    
        if verbose:
            
            
            if stdout:
                printBlue("STDOUT:\n"+stdout)
            if stderr:
                printBoldRed("STDERR:\n"+stderr)
        if not quiet:
            printGreen("Time taken:"+str(dt.timedelta(seconds=timeDiff)))
                
        exitCode=result.returncode
        
        ##Add to logs
        if logs:
            fullMessage=logMessage+"\n"+"exit code:"+str(exitCode)+"\texecution time:"+str(dt.timedelta(seconds=timeDiff))
            pyrpipeLoggerObject.cmdLogger.debug(fullMessage)
        
            ##log stdout
            pyrpipeLoggerObject.stdoutLogger.debug(logMessage+"\n"+stdout)
            ##log stderr
            pyrpipeLoggerObject.stderrLogger.debug(logMessage+"\n"+stderr)
            
            ##get the program used and log its path
            thisProgram=cmd[0]
            if thisProgram not in pyrpipeLoggerObject.loggedPrograms:
                ##get which thisProgram
                pyrpipeLoggerObject.envLogger.debug(thisProgram+":"+getProgramPath(thisProgram).strip())
                pyrpipeLoggerObject.loggedPrograms.append(thisProgram)
        
    
        if exitCode==0:
            return True
        return False
    
    except OSError as e:
        printBoldRed("Fatal Error occured"+str(e))
        return False
    except subprocess.CalledProcessError as e:
        printBoldRed("Fatal  Error occured:"+str(e))
        return False
    except:
        printBoldRed("Fatal Error occured during execution")
        return False
    



#modyfied from https://www.biostars.org/p/139422/
def isPairedSRA(pathToSraFile):
    """Function to test wheather a .sra file is paired or single.
    
    Parameters
    ----------
    arg1: string
        the path ro sra file
    """
    if not checkFilesExists(pathToSraFile):
        raise Exception("Error checking layout. {0} doesn't exist".format(pathToSraFile));
    
    try:
        fastqdCmd=["fastq-dump","-X","1","-Z","--split-spot", pathToSraFile]
        output = subprocess.check_output(fastqdCmd,stderr=subprocess.DEVNULL);
        numLines=output.decode("utf-8").count("\n")
        if(numLines == 4):
            return False;
        elif(numLines == 8):
            return True
        else:
            raise Exception("Unexpected output from fast-dump");
    except subprocess.CalledProcessError as e:
        raise Exception("Error running fastq-dump");
    
def getProgramPath(programName):
    whichCmd=['which',programName]
    out = subprocess.check_output(whichCmd,universal_newlines=True)
    return out
    
def checkDep(depList):
    """Check whether specified programs exist in the environment.
    This uses the which command to test whether a program is present.
    
    Parameters
    ----------
    arg1: list
        list of programs to test
        
    Returns
    -------
        bool True is all dependencies are satified, False otherwise.
    """
    errorFlag=False
    for s in depList:
        #printBlue("Checking "+s+"...")
        thisCmd=['which',s]
        if(getCommandReturnStatus(thisCmd)):
            #printGreen ("Found "+s)
            pass
        else:
            printBoldRed ("Can not find "+s)
            errorFlag=True
    if errorFlag:
        return False
    return True


#TODO: Re-implement following using native python libraries and move to utils
def deleteFileFromDisk(filePath):
    if checkFilesExists(filePath):
        rm_Cmd=['rm',filePath]
        rv= getCommandReturnStatus(rm_Cmd)
        return rv
    #if file doesn't exist return true
    return True

def deleteMultipleFilesFromDisk(*args):
    errorFlag=False
    for filePath in args:
        status=deleteFileFromDisk(filePath)
        if not status:
            errorFlag=True
    
    return not(errorFlag)

def moveFile(source,destination):
    """
    perform mv command
    """
    print("MOV:"+source+"-->"+destination)
    mv_cmd=['mv',source,destination]
    if not getCommandReturnStatus(mv_cmd):
        return False
    return True

#function to search files using find and return results as a list
#use os.scandir
def findFiles(path,name,recursive):
    if recursive:
        find_cmd=['find', path,'-type', 'f','-name',name]   
    else:
        find_cmd=['find', path,'-type', 'f','-maxdepth', '1','-name',name] 
    #print ("Executing: "+ ' '.join(find_cmd))
    #get output as string
    out = subprocess.check_output(find_cmd,universal_newlines=True)
    results=out.split()
    return results


if __name__ == "__main__": 
    print ("Logger")


