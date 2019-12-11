#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 11:16:40 2019

@author: usingh
"""

import sys
import os
import argparse
from pyrpipe import pyrpipe_utils as pu







def getCommandsFromLog(inFile):
    with open(inFile) as f:
        data=f.read().splitlines()
    commands=[]
    for l in data:
        if l.startswith("$"):
            commands.append(l.split("$ ")[1])
    return commands


"""
Valid options
help: print help

report: compile a full report of the analysis
shell: compile a shell script of all executed commands
failed: compile a list of failed commands
passed: compile a list of passed commands
benchmark: generate benchmarks from logs
"""

parser = argparse.ArgumentParser(prog='pyrpipe_diagnostic',add_help=True)
parser.add_argument("--verbose", help="increase output verbosity",action="store_true")
parser.add_argument("--report", help="increase output verbosity",action="store_true")
parser.add_argument("--shell", help="increase output verbosity",action="store_true")
parser.add_argument("--failed", help="increase output verbosity",action="store_true")
parser.add_argument("--passed", help="increase output verbosity",action="store_true")
parser.add_argument("--benchmark", help="increase output verbosity",action="store_true")
parser.add_argument("--logfile", help="increase output verbosity",action="store")


args = parser.parse_args()

logFile=args.logfile

if logFile is None:
    print("--logfile is required input. Exiting.")
    sys.exit(1)

print("logfile is:"+logFile)

#get log file basename
logFileDir=pu.getFileDirectory(logFile)
logFilename=pu.getFileName(logFile)
basename=logFilename.split(".")[0]

#print("Dir:{} fn:{} bn:{}".format(logFileDir,logFilename,basename))

#check all logs exist
cmdLog=os.path.join(logFileDir,basename+".CMD.log")
outLog=os.path.join(logFileDir,basename+".OUT.log")
errLog=os.path.join(logFileDir,basename+".ERR.log")
envLog=os.path.join(logFileDir,basename+".ENV.log")
if not pu.checkFilesExists(cmdLog,outLog,errLog,envLog):
    print("Please check missing log files. Exiting.")
    sys.exit(1)
    

###start processing
vFlag=args.verbose

if vFlag:
    print("pyrpipe diagnostic utility")


if args.report:
    print("Generating report")
    
    
if args.shell:
    print("Generating shell script")
    #read commands log
    commands=getCommandsFromLog(cmdLog)
    shebang="#!/bin/bash "
    #write to file
    outFile=basename+"_shell.sh"
    f=open(outFile,"w")
    f.write(shebang+"\n")
    f.write("\n".join(commands))
    
    if vFlag:
        print("\n\n".join(commands))
    print("shell commands written to:"+outFile)
    
    
    
    
    
    
    
    
    
    
    
    
    
    

