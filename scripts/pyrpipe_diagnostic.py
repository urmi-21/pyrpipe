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
import datetime as dt
#from markdownify import markdownify as md

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

from pyrpipe import report_templates


def parseEnvLog(envLog):
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
    return sysInfo,progList

def generateEnvReportTable(sysInfo,progList):
    """create html table to list environment
    Parameters
    ----------
    sysInfo: dict
        system information from env log file
    progList: dict
        programs and their information from env log
    """
    
    
    
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

def generatReportSummaryTable(sysInfo,progList):
    pass
    
    
def generateHTMLReport(templateFile,cmdLog,envLog,coverage='f'):
    """Generates html report
    Parameters
    ----------
    templatefile: string
        path to a template file
    cmdlog: string
        path to the log file
    envlog: string
        path to the env log file
    coverage: string
        tpye of report: full, summary, fail, pass
    
    """
    #vars for generating summary
    startTime=""
    endTime=""
    numCommands=0
    failedCommands=0
    passedCommands=0
    numPrograms=0
    progNames=[]
    #parse envLog
    sysInfo,progList=parseEnvLog(envLog)
    
    #get starttime #end time is calculated from log below
    startTime=dt.datetime.strptime(sysInfo['now'],"%y-%m-%d %H:%M:%S")
    #total progs used
    progNames=progList.keys()
    numPrograms=len(progNames)
    
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
    headHTML=pkg_resources.read_text(report_templates, 'head.html')
    #add file name
    headHTML+='\n<h2> <em>pyrpipe</em> report</h2>'
    headHTML+='\n<stron>file name:{}</strong>'.format(cmdLog)
    headHTML+='\n<hr><br><br>\n'
    
    fullHTML="\n<h2> Details </h2>"
    failColor="rgb(208,28,139)"
    passColor="rgb(77,172,38)" 
    for l in data:
        if not l.startswith("#"):
            thisDict=json.loads(l)
            numCommands+=1
            #add color to table
            if int(thisDict['exitcode'])==0:
                thisDict['statuscolor']=passColor
                passedCommands+=1
            else:
                thisDict['statuscolor']=failColor
                failedCommands+=1
            
            #program name
            programname=thisDict['cmd'].split(" ")[0]
            #add program version info
            newDict={**thisDict,**progList[programname]}
            
            #skip passed
            if coverage=='i' and int(thisDict['exitcode'])==0:
                continue
            #skip failed
            if coverage=='p' and int(thisDict['exitcode'])!=0:
                continue
            
            #escape all special html charecters
            for k, v in newDict.items():
                newDict[k] = escape(str(v))
            fullHTML=fullHTML+"\n"+template.render(newDict)
            
    #get start and runtime of last command
    lastDict=json.loads(data[-1])
    lastST=dt.datetime.strptime(lastDict['starttime'],"%y-%m-%d %H:%M:%S")
    lastruntime= dt.datetime.strptime(lastDict['runtime'],"%H:%M:%S")
    deltaTime = dt.timedelta(hours=lastruntime.hour, minutes=lastruntime.minute, seconds=lastruntime.second)
    endTime=lastST+deltaTime
    
    
            
    #generate summary
    summary='\n<h2> Summary </h2>'
    summary+='\n<div class="summary">'
    summary+='\n<code>Time start: {}     Time end: {}      Total time: {}</code>'.format(str(startTime),str(endTime), str(endTime-startTime))
    summary+='\n<br><br>'
    summary+='\n<span>Num commands: {}</span>'.format(numCommands)
    summary+='\n<br><br>'
    summary+='\nNum failed commands: {}'.format(failedCommands)
    summary+='\n<br><br>'
    summary+='\nNum passed commands: {}'.format(passedCommands)
    summary+='\n<br><br>'
    summary+='\nTotal programs: {}'.format(numPrograms)
    summary+='\n<br><br>'
    summary+='\nNum passed commands: {}'.format(",".join(progNames))
    summary+='\n<br><br>'
    summary+='\n</div>'
        
    envTable=generateEnvReportTable(sysInfo,progList)
    headHTML=headHTML+summary+fullHTML+envTable+"\n</body>\n</html>"
    return headHTML
    

def writeHtmlToPdf(htmlText,outFile):
    if not outFile.endswith(".pdf"):
        outFile=outFile+".pdf"
   
    #read css
    cssFile = pkg_resources.read_text(report_templates, 'simple.css')
    #HTML(string=htmlText).write_pdf(outFile)
    #HTML(string=htmlText).write_pdf(outFile, stylesheets=[CSS('/home/usingh/work/urmi/hoap/pyrpipe/pyrpipe/report_templates/simple.css')])
    HTML(string=htmlText).write_pdf(outFile, stylesheets=[CSS(string=cssFile)])
    print("Report written to {}".format(outFile))
    
def writeHtml(htmlText,outFile):
    if not outFile.endswith(".html"):
        outFile=outFile+".html"
    #read css
    #cssFile = pkg_resources.read_text(report_templates, 'simple.css')
    #HTML(string=htmlText).write_pdf(outFile)
    #HTML(string=htmlText).write_pdf(outFile, stylesheets=[CSS('/home/usingh/work/urmi/hoap/pyrpipe/pyrpipe/report_templates/simple.css')])
    #HTML(string=htmlText).write_pdf(outFile, stylesheets=[CSS(string=cssFile)])
    f=open(outFile,'w')
    f.write(htmlText)
    f.close()
    
    print("Report written to {}".format(outFile))
    
#TODO: write markdown report
def writeHtmlToMarkdown(htmlText,outFile):
    pass
    #if not outFile.endswith(".md"):
    #    outFile=outFile+".md"
    #print(md(htmlText))

def getCommandsFromLog(inFile):
    with open(inFile) as f:
        data=f.read().splitlines()
    commands=[]
    for l in data:
        if not l.startswith("#"):
            commands.append(json.loads(l)["cmd"])
    return commands

def generateBashScript(logFile,outFile,filterList,coverage='a'):
    commands=getCommandsFromLog(logFile,filterList,coverage)
    shebang="#!/bin/bash "
    #write to file
    outFile=basename+"_shell.sh"
    f=open(outFile,"w")
    f.write(shebang+"\n")
    f.write("\n".join(commands))
    
    if vFlag:
        print("\n\n".join(commands))
    print("shell commands written to:"+outFile)

"""
Valid options
help: print help

report: compile a full report of the analysis
bash: compile a shell script of all executed commands
failed: compile a list of failed commands
passed: compile a list of passed commands
benchmark: generate benchmarks from logs
"""

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
"""

def checkEnvLog(logFile):
    #check all logs exist
    logFileDir=pu.getFileDirectory(logFile)
    basename=pu.getFileBaseName(logFile)
    envLog=os.path.join(logFileDir,basename+"ENV.log")
    if not pu.checkFilesExists(logFile,envLog):
        print("Please check missing log files. Exiting.")
        sys.exit(1)
    return envLog


def report():
    
    parser = argparse.ArgumentParser(
   
            description='pyrpipe diagnostic utility\nGenerate analysis report.',
            
            usage='''pyrpipe_diagnostic report [<args>] <logfile>
                    
                    ''')    
    parser.add_argument('-o', help='out file \ndefault: same as input logfile',action="store")
    parser.add_argument('-e', help='report output type: [md,pdf,html] \ndefault: pdf',default='pdf',action="store")
    parser.add_argument('-c',help='Report options [(f)ull,fa(i)l,(p)ass]\ndefault: f',default='f',action="store")
    parser.add_argument('-v',help='verbose',action="store_true")
    parser.add_argument('logfile', help='The log file generated by pyrpipe',action="store")
    args = parser.parse_args(sys.argv[2:])
    
    logFile=args.logfile
    envLog=checkEnvLog(logFile)    
    #parse args
    vFlag=args.v
    if vFlag:
        print("Generating report")
    outFile=""
    if args.o is None:
        outFile=pu.getFileBaseName(args.logfile)
    else:
        outFile=args.o
    outFile+='.'+args.e
    
    
    if args.e in ['pdf','html']:
        htmlReport=generateHTMLReport('simpleDiv.html',logFile,envLog,coverage=args.c)
        if args.e=='pdf':
            writeHtmlToPdf(htmlReport,outFile)
        else:
            writeHtml(htmlReport,outFile)
    elif args.e == 'md':
        pass
    else:
        pu.printBoldRed("unknown extension:"+args.e+". Exiting")
    
    
    
def bash():
    print("Generating bash script")
    parser = argparse.ArgumentParser(
   
            description='pyrpipe diagnostic utility\nGenerate bash script.',
            
            usage='''pyrpipe_diagnostic report [<args>] <logfile>
                    
                    ''')    
    parser.add_argument('-o', help='out file \ndefault: same as input logfile',action="store")
    parser.add_argument('-c',help='Dump command options [(a)ll,fa(i)l,(p)ass]\ndefault: a',default='a',action="store")
    parser.add_argument('-v',help='verbose',action="store_true")
    parser.add_argument('-f',help='Filter by programs. Provide a comma-separated list e.g., prefetch,STAR,bowtie2 \ndefault None')
    parser.add_argument('logfile', help='The log file generated by pyrpipe',action="store")
    args = parser.parse_args(sys.argv[2:])
    
    logFile=args.logfile  
    #parse args
    vFlag=args.v
    if vFlag:
        print("Generating report")
    outFile=""
    if args.o is None:
        outFile=pu.getFileBaseName(logFile)
    else:
        outFile=args.o
    outFile+='.sh'
    
    filters=[]
    if args.f is not None:
        filters= args.f.split(',')
    
    generateBashScript(logFile,outFile,filters,args.c)
    

def benchmark():
    print("Generating benchmarks")
    parser = argparse.ArgumentParser(
   
            description='pyrpipe diagnostic utility\nGenerate benchmark report.',
            
            usage='''pyrpipe_diagnostic report [<args>] <logfile>
                    
                    ''')    
    parser.add_argument('-o', help='out file \ndefault: same as input logfile',action="store")
    parser.add_argument('-e', help='report output type: [MD,PDF,HTML] \ndefault: PDF',default='PDF',action="store")
    parser.add_argument('-v',help='verbose',action="store_true")
    parser.add_argument('-f',help='Filter by programs. Provide a comma-separated list e.g., prefetch,STAR,bowtie2 \ndefault None')
    parser.add_argument('logfile', help='The log file generated by pyrpipe',action="store")
    args = parser.parse_args(sys.argv[2:])
    print(str(args))



subcommands=['report','bash','benchmark','all']
parser = argparse.ArgumentParser(
            
            description='pyrpipe diagnostic utility',
            
            usage='''pyrpipe_diagnostic <command> [<args>] <logfile>
                    The commands are:
                    report     Generate analysis report
                    bash      Generate all commands to bash script
                    benchmark Generate bemchmarks
                    ''')
parser.add_argument('command', help='Subcommand to run [report,bash,benchmark,all]')


#parse first and last argument as subcommand and logfile 
args = parser.parse_args(sys.argv[1:2])
if args.command not in subcommands:
    print ('Unrecognized command')
    parser.print_help()
    exit(1)

if args.command == 'report':
    report()
elif args.command == 'bash':
    bash()
elif args.command == 'benchmark':
    benchmark()





exit(1)
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    

