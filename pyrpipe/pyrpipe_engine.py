#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 16:26:02 2019

@author: usingh

Methods and classes related to execution and logging
"""

import os
import subprocess
import datetime as dt
import time

#for log
import time
from datetime import datetime 
from datetime import timedelta
import logging
import os
#from sinfo import sinfo
import platform
from multiprocessing import cpu_count

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

class pyrpipeLogger():
    def __init__(self,name,logfile,formatter,level=logging.DEBUG):
        self.__name__="pyrpipeLogger"
        self.logger=self.createLogger(name,logfile,formatter,level=logging.DEBUG)
    
    def createLogger(name,logfile,formatter,level=logging.DEBUG):
        #Get different loggers
        handler = logging.FileHandler(logfile)        
        handler.setFormatter(formatter)
        
        logger = logging.getLogger(name)
        logger.setLevel(level)
        logger.addHandler(handler)
        return logger
        
###create loggers
        
timestamp=str(datetime.now()).split(".")[0].replace(" ","-").replace(":","_")

#create logs directory

logsDir=os.path.join(os.getcwd(),"pyrpipe_logs")
if not os.path.isdir(logsDir):
    os.mkdir(logsDir)
#log names
cmdLogFname=os.path.join(logsDir,timestamp+"_pyrpipeCMD.log")
cmdLogOb=pyrpipeLogger("cmd",cmdLogFname,LogFormatter())
cmdLogOb.logger.debug("#START LOGNEWWWWW")
    
    
    
    
    
    
    
    
    
##############Functions###########################
    

def getCommandReturnValue(cmd):
    result = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    stdout,stderr = result.communicate()
    return result.returncode

def getCommandReturnStatus(cmd):
    returnValue=getCommandReturnValue(cmd)
    if returnValue==0:
        return True
    return False

#prints stdout in real time. optimal for huge stdout and no stderr
def executeCommandOld(cmd):
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


def executeCommand(cmd,verbose=False,quiet=False):
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
        fullMessage=logMessage+"\n"+"exit code:"+str(exitCode)+"\texecution time:"+str(dt.timedelta(seconds=timeDiff))
        pl.commandLogger.debug(fullMessage)
        
        ##log stdout
        pl.stdOutLogger.debug(logMessage+"\n"+stdout)
        ##log stderr
        pl.stdErrLogger.debug(logMessage+"\n"+stderr)
        
        ##get the program used and log its path
        thisProgram=cmd[0]
        if thisProgram not in pl.loggedPrograms:
            ##get which thisProgram
            pl.envLogger.debug(thisProgram+":"+getProgramPath(thisProgram).strip())
            pl.loggedPrograms.append(thisProgram)
            
    
        if exitCode==0:
            return True
        return False
    
    except OSError as e:
        print("Fatal Error occured"+str(e))
        return False
    except subprocess.CalledProcessError as e:
        print("Fatal  Error occured:"+str(e))
        return False
    except:
        print("Fatal Error occured")
        return False
    
    
    
    
#function to search files using find and return results as a list
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
        if(getCommandReturnValue(thisCmd)==0):
            #printGreen ("Found "+s)
            pass
        else:
            printBoldRed ("Can not find "+s)
            errorFlag=True
    if errorFlag:
        return False
    return True


def deleteFileFromDisk(filePath):
    if checkFilesExists(filePath):
        rm_Cmd=['rm',filePath]
        rv= getCommandReturnStatus(rm_Cmd)
        return rv
    #if file doesn't exist return true
    return True


def moveFile(source,destination):
    """
    perform mv command
    """
    print("MOV:"+source+"-->"+destination)
    mv_cmd=['mv',source,destination]
    if not getCommandReturnStatus(mv_cmd):
        return False
    return True
