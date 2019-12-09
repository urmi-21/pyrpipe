#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 12:04:28 2019

@author: usingh
"""

import os
import subprocess
import dill
import datetime as dt
import time
from pyrpipe import pyrpipe_logger as pl



#functions to print in color

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def printBoldRed(text):
    print (bcolors.FAIL + bcolors.BOLD+ text + bcolors.ENDC)

def printGreen(text):
    print (bcolors.OKGREEN + text + bcolors.ENDC)

def printBlue(text):
    print (bcolors.OKBLUE + text + bcolors.ENDC) 

######End color functions###################

def getTimestamp(shorten=False):
    
    timestamp=str(dt.datetime.now()).split(".")[0].replace(" ","-")
    if shorten:
        timestamp=timestamp.replace("-","").replace(" ","").replace(":","")
    return timestamp
    

def savePyrpipeWorkspace(filename="myWorkspace",outDir=""):
    """Save current workspace using dill.
    """
    #timestamp format YYYYMMDDHHMISE
    timestamp=getTimestamp(True)
    
    
    if not outDir:        
        outDir=os.getcwd()
    
    outFile=os.path.join(outDir,filename)
    outFile=outFile+"_"+timestamp+".pyrpipe"
    
    #save workspace
    dill.dump_session(outFile)
    print("Session saved.")


def restorePyrpipeWorkspace(file):
    if not checkFilesExists(file):
        print(file+" doesn't exist")
        return False
    #load the session
    dill.load_session(file)
    print("Session restored.")
    return True

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


def executeCommand(cmd,verbose=False):
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
            printBlue(logMessage)
            
            if stdout:
                printBlue("STDOUT:\n"+stdout)
            if stderr:
                printBoldRed("STDERR:\n"+stderr)
    
            printBlue("Time taken:"+str(dt.timedelta(seconds=timeDiff)))
                
        exitCode=result.returncode
        
        ##Add to logs
        fullMessage=logMessage+"\n"+"exit code:"+str(exitCode)+"\texecution time:"+str(dt.timedelta(seconds=timeDiff))
        pl.commandLogger.debug(fullMessage)
        
        ##get the program used
    
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
        print("Fatal Error occured:"+str(e))
        return False


def getSRADownloadPath(srrID):
    if len(srrID) <6:
        return None
    
    parentPath="anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/"
    parentPath=parentPath+srrID[0:3]+"/"+srrID[0:6]+"/"+srrID+"/"+srrID+".sra"
    
    return parentPath
    
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

def checkPathsExists(*args):
    """
    check if a directory exists
    """
    failFlag=False
    for path in args:
        if not (os.path.exists(path) and os.path.isdir(path)):
            failFlag=True
    if failFlag==True:
        return False
    return True


def checkFilesExists(*args):
    failFlag=False
    for path in args:
        if not os.path.isfile(path):
            #printBoldRed("File not found: "+path)
            failFlag=True
    
    if failFlag:
        return False
    return True

def checkHisatIndex(index):
    return checkFilesExists(index+".1.ht2")

def checkBowtie2Index(index):
    return checkFilesExists(index+".1.bt2")
    

def bytetoReadable(sizeInBytes):
    """
    function to convert bytes to human readable format (MB,GB ...)
    """
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if sizeInBytes < 1024.0:
            return "%3.1f %s" % (sizeInBytes, x)
        sizeInBytes /= 1024.0


def getFileSize(file_path):
    """
    Return File size in human readable format
    """
    
    if (checkFilesExists(file_path)):
        file_info = os.stat(file_path)
        return bytetoReadable(file_info.st_size)
    

#TODO: override in case of empty list
def parseJavaStyleArgs(validArgsList,passedArgs):
    """
    Function creates arguments to pass to unix systems through popen
    Parameters
    ----------
    arg1 : list
        list of valid arguments. Invalid arguments will be ignored
    arg2: keyword value argument list to be parsed
        
    Returns
    -------
        list
            a list with command line arguments to be used with popen

        Examples
        --------
        >>> parseJavaStyleArgs(['A','B','-C'], {"A": "3", "B": "22","-C":""})
        ['A=3', 'B=22', '-C']
    """
    popenArgs=[]
    specialArgList=["--"]
    appendAtEndArgs=[]
    
    #empty list supplied consider all armunets valid
    if len(validArgsList)<1:
        validArgsList=passedArgs.keys()
        #above command will also add specialArgs, remove those
        for x in specialArgList:
            if x in validArgsList:
                validArgsList.remove(x)
    
    for key, value in passedArgs.items():
        #check if key is a valid argument
        if key in validArgsList:
            #do not add emty parameters e.g. -q or -v
            if len(value)>0:
                popenArgs.append(key+"="+value)
            else:
                popenArgs.append(key)
        elif key in specialArgList:
            appendAtEndArgs.extend(value)
        
        else:
            print("Unknown argument {0} {1}. ignoring...".format(key, value))
    popenArgs.extend(appendAtEndArgs)
    return popenArgs
    
    
def parseUnixStyleArgs(validArgsList,passedArgs):
    """
    Function creates arguments to pass to unix systems through popen
    Parameters
    ----------
    arg1 : list
        list of valid arguments. Invalid arguments will be ignored
    arg2: keyword value argument list to be parsed
        
    Returns
    -------
        list
            a list with command line arguments to be used with popen

        Examples
        --------
        >>> parseUnixStyleArgs(['-O','-t','-q'], {"-O": "./test", "Attr2": "XX","--":("IN1","IN2")})
        Unknown argument Attr2 XX. ignoring...
        ['-O', './test', 'IN1', 'IN2']
    """
    
        
    
    popenArgs=[]
    """
    Define some special arguments.
    -- to pass input which don't follow any flag e.g. mypro.sh input1 input2
    """
    specialArgList=["--"]
    appendAtEndArgs=[]
    
    #empty list supplied consider all armunets valid
    if len(validArgsList)<1:
        validArgsList=passedArgs.keys()
        #above command will also add specialArgs, remove those
        for x in specialArgList:
            if x in validArgsList:
                validArgsList.remove(x)
        
        
    for key, value in passedArgs.items():
        #check if key is a valid argument
        if key in validArgsList:
            popenArgs.append(key)
            #do not add emty parameters e.g. -q or -v
            if len(value)>0:
                    popenArgs.append(value)
        elif key in specialArgList:
            appendAtEndArgs.extend(value)
        else:
            print("Unknown argument {0} {1}. ignoring...".format(key, value))
    popenArgs.extend(appendAtEndArgs)
    return popenArgs
    

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
        printBlue("Checking "+s+"...")
        thisCmd=['which',s]
        if(getCommandReturnValue(thisCmd)==0):
            printGreen ("Found "+s)
        else:
            printBoldRed ("Can not find "+s)
            errorFlag=True
    if errorFlag:
        return False
    return True


def getFileDirectory(filePath):
    return os.path.split(filePath)[0]

def getFileName(filePath):
    return os.path.split(filePath)[1]

def getFileExtension(filePath):
    return os.path.splitext(filePath)[1]

def getFileBaseName(filePath):
    """return file name without the extension
    arg1:
        file path or file name
    """
    return os.path.splitext(getFileName(filePath))[0]

def deleteMultipleFilesFromDisk(*args):
    errorFlag=False
    for filePath in args:
        status=deleteFileFromDisk(filePath)
        if not status:
            errorFlag=True
    
    return not(errorFlag)
    
def deleteFileFromDisk(filePath):
    if checkFilesExists(filePath):
        rm_Cmd=['rm',filePath]
        rv= getCommandReturnStatus(rm_Cmd)
        return rv
    #if file doesn't exist return true
    return True


def mkdir(dirPath):
    """Create a directory
    """
    print("creating Dir:"+dirPath)
    try:
        os.mkdir(dirPath)
    except OSError:
        return False
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
    

if __name__ == "__main__":
    print("main")
    #print(parseUnixStyleArgs(['-O','-t','-q'], {"-O": "./test", "Attr2": "XX","--":("IN1","IN2")}))
    print(parseJavaStyleArgs(['A','B','-C'], {"A": "3", "B": "22","-C":"","--":("-Xmx2g","-da",)}))
