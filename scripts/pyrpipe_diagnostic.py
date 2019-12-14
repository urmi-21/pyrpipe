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
from html import escape
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

def generateEnvReportTable(sysInfo,progList):
    """create html table to list environment
    """
    i=0
    os=sysInfo['os']
    now=sysInfo['now']
    sysmodulesList=sysInfo['sysmodules'].strip('][').split(', ')
    syspathList=sysInfo['syspath'].strip('][').split(', ')
    python=sysInfo['python']
    cpu=sysInfo['cpu']
    
    ##create system info table
    tabStr='\n<h2>Environment Information</h2>'
    tabStr+='\n <div class="envtabs">'
    #add program Info
    tabStr+='\n<table class="programinfo" >'
    tabStr+='<tr><th colspan="3">Programs</th></tr>'
    tabStr+='\n<tr> <th>{}</th> <th>{}</th> <th>{}</th>  </tr>'.format("name","version","path")
    for k, v in progList.items():
        tabStr+='\n<tr> <td>{}</td> <td>{}</td> <td>{}</td>  </tr>'.format(v["name"],v["version"],v["path"])
    tabStr+='\n</table>'    
    tabStr+='\n<br><br>'
    tabStr+='\n<table class="sysInfotable" >'
    tabStr+='<tr><th colspan="2">System Information</th></tr>'
    tabStr+='\n<tr>    <td>Time at collection</td> <td>{}</td>    </tr>'.format(now)
    tabStr+='\n<tr>    <td>Python</td> <td>{}</td>    </tr>'.format(python)
    tabStr+='\n<tr>    <td>Operating system</td> <td>{}</td>    </tr>'.format(os)
    tabStr+='\n<tr>    <td>CPU</td> <td>{}</td>    </tr>'.format(cpu)
    tabStr+='\n</table>'
    tabStr+='\n<br><br>'
    #add sys modules table
    tabStr+='\n<table class="sysmodules" >'
    tabStr+='<tr><th colspan="1">sys.modules</th></tr>'
    for s in sysmodulesList:
        tabStr+='\n<tr><td>{}</td></tr>'.format(s)
    tabStr+='\n</table>'
    tabStr+='\n<br><br>'
    #add sys path table
    tabStr+='\n<table class="syspath" >'
    tabStr+='<tr><th colspan="1">sys.path</th></tr>'
    for s in syspathList:
        tabStr+='\n<tr><td>{}</td></tr>'.format(s)
    tabStr+='\n</table>'
    tabStr+='\n<br><br>'
    tabStr+='\n</div>'
            
    
    
    return tabStr
    
    
def generateHTMLReport(templateFile,cmdLog,envLog):
    
    #parse the env log
    with open(envLog) as f:
        envdata=f.read().splitlines()
    
    sysInfo={}
    progList={}
    for l in envdata:
        if not l.startswith("#"):
            if not sysInfo:
                sysInfo=json.loads(l)
            else:
                thisProgram=json.loads(l)
                progList[thisProgram['name']]=thisProgram
                
                
    #print("SYSINFO:"+str(sysInfo))
    #print("PROGLIST:"+str(progList))
    
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
    #add file name
    fullHTML+='\n<h2> <em>pyrpipe</em> report</h2>'
    fullHTML+='\n<stron>file name:{}</strong>'.format(cmdLog)
    fullHTML+='\n<hr><br><br>\n'
    
    failColor="rgb(208,28,139)"
    passColor="rgb(77,172,38)" 
    for l in data:
        if not l.startswith("#"):
            thisDict=json.loads(l)
            #add color to table
            if int(thisDict['exitcode'])==0:
                thisDict['statuscolor']=passColor
            else:
                thisDict['statuscolor']=failColor
            
            #program name
            programname=thisDict['cmd'].split(" ")[0]
            #add program version info
            newDict={**thisDict,**progList[programname]}
            #escape all special html charecters
            for k, v in newDict.items():
                newDict[k] = escape(str(v))
            
            fullHTML=fullHTML+"\n"+template.render(newDict)#return html
            
            
    envTable=generateEnvReportTable(sysInfo,progList)
    fullHTML=fullHTML+envTable+"\n</body>\n</html>"
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
        print(htmlReport)
        pass
        
    
    writeHtmlToPdf(htmlReport,logFilename)
    #writeHtmlToMarkdown(htmlReport,"rep")
    
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    

