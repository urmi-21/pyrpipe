#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 12:04:28 2019

@author: usingh
"""

import os
import subprocess

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

def getCommandReturnValue(cmd):
	result = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
	stdout,stderr = result.communicate()
	return result.returncode

def executeCommand(cmd):
      
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
            yield stdout_line 
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)
        

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
    print ("Executing: "+ ' '.join(find_cmd))
    #get output as string
    out = subprocess.check_output(find_cmd,universal_newlines=True)
    results=out.split()
    return results

def checkPathsExists(*args):
	failFlag=False
	for path in args:
		if not os.path.exists(path):
			printBoldRed("Path not found: "+path)
			failFlag=True
	if failFlag==True:
		return 1
	return 0


def checkFilesExists(*args):
	failFlag=False
	for path in args:
		if not os.path.isfile(path):
			printBoldRed("File not found: "+path)
			failFlag=True
	
	if failFlag:
		return 1
	return 0

def checkHisatIndex(index):
	return checkFilesExists(index+".1.ht2")
	
    

if __name__ == "__main__":
    #test
    #print(getSRADownloadPath('SRR002328'))
    print(findFiles("/home/usingh/work/urmi","*.py",False))
 

    
