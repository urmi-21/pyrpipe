#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 11:16:40 2019

@author: usingh
"""

import sys
import os
import argparse
import json
from pyrpipe import pyrpipe_utils as pu
from jinja2 import Environment, BaseLoader
from weasyprint import HTML
try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

from pyrpipe import report_templates




def getCommandsFromLog(inFile):
    with open(inFile) as f:
        data=f.read().splitlines()
    commands=[]
    for l in data:
        if not l.startswith("#"):
            commands.append(json.loads(l)["cmd"])
    return commands

def generateHTMLReport(templateFile,cmdLog,envLog):
    
    #read the template
    templateHTMLFile = pkg_resources.read_text(report_templates, templateFile)
    #load template file
    template = Environment(loader=BaseLoader()).from_string(templateHTMLFile)
    
    """
    for each command executed get the report parameters.
    Template should specify the parameters:
        cmd: the command executed
        stdout: the stdout
        stderr: the stderr
        exitcode: return code obtained
        starttime: time started
        runtime: execution time
        
        example record:{'cmd':logMessage,
                 'exitcode':str(exitCode),
                 'runtime':str(dt.timedelta(seconds=timeDiff)),
                 'starttime':str(strStartTime),
                 'stdout':stdout,
                 'stderr':stderr                 
                }
    """
    with open(cmdLog) as f:
        data=f.read().splitlines()
    #save log for each command
    logs=[]
    fullHTML=""
    for l in data:
        if not l.startswith("#"):
            thisDict=json.loads(l)
            fullHTML=fullHTML+template.render(thisDict)#return html
    
    return fullHTML
    
        


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
parser.add_argument("--log", help="increase output verbosity",action="store")


args = parser.parse_args()

logFile=args.log

if logFile is None:
    print("--logfile is required input. Exiting.")
    sys.exit(1)

print("logfile is:"+logFile)

#get log file basename
logFileDir=pu.getFileDirectory(logFile)
logFilename=pu.getFileName(logFile)
basename=pu.getFileBaseName(logFile)

#print("Dir:{} fn:{} bn:{}".format(logFileDir,logFilename,basename))

#check all logs exist
envLog=os.path.join(logFileDir,basename+"ENV.log")
if not pu.checkFilesExists(logFile,envLog):
    print("Please check missing log files. Exiting.")
    sys.exit(1)
    

###start processing
vFlag=args.verbose

if vFlag:
    print("pyrpipe diagnostic utility")


if args.report:
    print("Generating report")
    
    #env = Environment(loader=FileSystemLoader('.'))
    """
    temp = pkg_resources.read_text(report_templates, 'simple.html')
    template = Environment(loader=BaseLoader()).from_string(temp)
    template_vars = {"title" : "Sales Funnel Report - National",
                 "national_pivot_table": "RANSJKDAJKDSA"}
    html_out = template.render(template_vars)
    print(html_out)
    HTML(string=html_out).write_pdf("report.pdf")
    
    temp = pkg_resources.read_text(report_templates, 'block.html')
    template = Environment(loader=BaseLoader()).from_string(temp)
    
    fullH=""
    
    for i in range(1,10):
        template_vars = {"title" : str(i)}
        fullH=fullH+template.render(template_vars)
    
    print (fullH)
    HTML(string=fullH).write_pdf("report2.pdf")
    """
    
    htmlReport=generateHTMLReport('basic.html',logFile,envLog)
    
    print(htmlReport)
    
    
    
if args.shell:
    print("Generating shell script")
    #read commands log
    commands=getCommandsFromLog(logFile)
    shebang="#!/bin/bash "
    #write to file
    outFile=basename+"_shell.sh"
    f=open(outFile,"w")
    f.write(shebang+"\n")
    f.write("\n".join(commands))
    
    if vFlag:
        print("\n\n".join(commands))
    print("shell commands written to:"+outFile)
    
    
    
    
    
    
    
    
    
    
    
    
    
    

