#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 12:04:28 2019

@author: usingh
"""

import os
import datetime as dt



#functions to print in color

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    LightMagenta = "\033[95m"
    LightYellow  = "\033[93m"
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

def printInfoMessage(text):
    print (bcolors.LightMagenta + text + bcolors.ENDC) 
def printYellow(text):
    print (bcolors.LightYellow + text + bcolors.ENDC) 

######End color functions###################

def getTimestamp(shorten=False):
    
    timestamp=str(dt.datetime.now()).split(".")[0].replace(" ","-")
    if shorten:
        timestamp=timestamp.replace("-","").replace(" ","").replace(":","")
    return timestamp
    
"""
moved to session.py
def savePyrpipeWorkspace(filename="myWorkspace",outDir=""):
    #Save current workspace using dill.
    
    #timestamp format YYYYMMDDHHMISE
    timestamp=getTimestamp(True)
    
    
    if not outDir:        
        outDir=os.getcwd()
    
    outFile=os.path.join(outDir,filename)
    outFile=outFile+"_"+timestamp+".pyrpipe"
    
    #save workspace
    dill.dump_session(outFile)
    print("Session saved to: "+outFile)


def restorePyrpipeWorkspace(file):
    if not checkFilesExists(file):
        print(file+" doesn't exist")
        return False
    #load the session
    dill.load_session(file)
    print("Session restored.")
    return True
"""




def getSRADownloadPath(srrID):
    if len(srrID) <6:
        return None
    
    parentPath="anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/"
    parentPath=parentPath+srrID[0:3]+"/"+srrID[0:6]+"/"+srrID+"/"+srrID+".sra"
    
    return parentPath
    


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
    



def mkdir(dirPath):
    """Create a directory
    """
    print("creating Dir:"+dirPath)
    try:
        os.mkdir(dirPath)
    except OSError:
        return False
    return True


    

if __name__ == "__main__":
    print("main")
    #print(parseUnixStyleArgs(['-O','-t','-q'], {"-O": "./test", "Attr2": "XX","--":("IN1","IN2")}))
    print(parseJavaStyleArgs(['A','B','-C'], {"A": "3", "B": "22","-C":"","--":("-Xmx2g","-da",)}))
