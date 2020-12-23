#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 13:23:55 2020

@author: usingh
"""

from pyrpipe import pyrpipe_engine as pe
from pyrpipe import reports

logfile='tests/test_files/pyrpipe_logs/2020-01-22-18_14_47_pyrpipe.log'
envlogfile='tests/test_files/pyrpipe_logs/2020-01-22-18_14_47_pyrpipeENV.log'
outdir='tests/testout/testreport'
logdir='tests/test_files/pyrpipe_logs'

def test_report():
    cmd="pyrpipe_diagnostic report -o tests/testout/testreport tests/test_files/pyrpipe_logs/2020-01-22-18_14_47_pyrpipe.log"
    st=pe.execute_command(cmd.split(),objectid="",command_name="")
    assert st==True, "report failed"


def test_benchmark():
    cmd="pyrpipe_diagnostic benchmark -t tests/testout/bmeport tests/test_files/pyrpipe_logs/2020-01-22-18_14_47_pyrpipe.log"
    st=pe.execute_command(cmd.split(),objectid="",command_name="")
    assert st==True, "benchmark failed"
    
def test_shell():
    cmd="pyrpipe_diagnostic shell -o tests/testout/shellcmds tests/test_files/pyrpipe_logs/2020-01-22-18_14_47_pyrpipe.log"
    st=pe.execute_command(cmd.split(),objectid="",command_name="")
    assert st==True, "shell failed"
    
    
def test_multiqc():
    cmd="pyrpipe_diagnostic multiqc -o tests/testout/mqcreport -t tests/testout/mqctmp tests/test_files/pyrpipe_logs/2020-01-22-18_14_47_pyrpipe.log"
    st=pe.execute_command(cmd.split(),objectid="",command_name="")
    assert st==True, "multiqc failed"
    
def test_multiqc2():
    cmd="pyrpipe_diagnostic multiqc -o tests/testout/mqcreport -t tests/testout/mqctmp tests/test_files/pyrpipe_logs"
    st=pe.execute_command(cmd.split(),objectid="",command_name="")
    assert st==True, "multiqc2 failed"
    
def test_multiqc3():    
    reports.generate_multiqc(directory='tests/testout',tempDir='tests/testout/mqctmp3',outDir='tests/testout/mqcout3')
      
    
def test_report2():
    reports.generateBenchmarkReport(logfile,envlogfile,None,'tests/testout/mqctmp',outFile="",verbose=False)

def test_summary():
    reports.generate_summary(logfile,envlogfile)
    

def test_checklog():
    assert reports.checkEnvLog(logfile) !=None,'check env log fail'
    
def test_generate_bash():
    reports.generateBashScript(logfile,logfile+'sh',None,coverage='a')
    reports.generateBashScript(logfile,logfile+'sh',None,coverage='i')
    reports.generateBashScript(logfile,logfile+'sh',None,coverage='p')
    
def test_get_stdout():
    assert reports.getStdoutFromLog(logfile,None,'a') != None, 'get stdout fail'
    
def test_HTMLReport():
    reports.generateHTMLReport('simpleDiv.html',logfile,envlogfile,coverage='f')
    
    
def test_multiqc4():
    cmd="rm -r tests/testout"
    assert pe.execute_command(cmd.split(),objectid="",command_name="")==True,'removing failed'
    #reports.generate_multiqc_from_log(logfile,None,'tests/testout/mqctmp4',outDir='tests/testout/mqcout4')   
    
    
    
    