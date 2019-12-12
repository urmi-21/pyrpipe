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
from weasyprint import HTML,CSS
#from markdownify import markdownify as md

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
    
    #read head.html
    fullHTML=pkg_resources.read_text(report_templates, 'head.html')
    #print(fullHTML)
    
    for l in data:
        if not l.startswith("#"):
            thisDict=json.loads(l)
            fullHTML=fullHTML+"\n"+template.render(thisDict)#return html
            
            
    
    fullHTML=fullHTML+"\n</body>\n</html>"
    return fullHTML
    

def writeHtmlToPdf(htmlText,outFile):
    if not outFile.endswith(".pdf"):
        outFile=outFile+".pdf"
   
    #read css
    cssFile = pkg_resources.read_text(report_templates, 'simple.css')
    #HTML(string=htmlText).write_pdf(outFile)
    #HTML(string=htmlText).write_pdf(outFile, stylesheets=[CSS('/home/usingh/work/urmi/hoap/pyrpipe/pyrpipe/report_templates/simple.css')])
    HTML(string=htmlText).write_pdf(outFile, stylesheets=[CSS(string=cssFile)])
    print("Report written to {}".format(outFile))
    
#TODO: write markdown report
def writeHtmlToMarkdown(htmlText,outFile):
    pass
    #if not outFile.endswith(".md"):
    #    outFile=outFile+".md"
    #print(md(htmlText))


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
    
    htmlReport=generateHTMLReport('simpleDiv.html',logFile,envLog)
    
    if vFlag:
        pass
        #print(htmlReport)
    
    writeHtmlToPdf(htmlReport,"rep")
    writeHtmlToMarkdown(htmlReport,"rep")
    
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    

