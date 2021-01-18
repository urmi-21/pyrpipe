#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 16:26:02 2019

@author: usingh

Methods and classes related to execution and logging
"""


import subprocess
import time
from datetime import timedelta
import sys
from pyrpipe import pyrpipe_utils as pu
import json
from shutil import which

#import dryrun from Conf
from pyrpipe import _dryrun
from pyrpipe import _safe
from pyrpipe import _logs_dir
from pyrpipe import _log_name
from pyrpipe import _logging
from pyrpipe import _verbose
from pyrpipe import PyrpipeLogger

###create logger


if _logging:
    pyrpipe_logger=PyrpipeLogger(_log_name,logdir=_logs_dir)
    pu.print_yellow("Logs will be saved to {}".format(pyrpipe_logger.logger_path))
    
"""
All functions that interact with shell are defined here. 
"""
    
##############Functions###########################

#a decorator to skip functions in safe mode
def skippable(func):
    """
    Skip a function execution in safemode
    """
    if not _safe:
        return func
    
    return lambda *args,**kwargs: True
    

#decorator function for dry runs
def dryable(func):
    """
    decorator function for drying all functions capable of executing commands
    """
    if not _dryrun:
        return func
    def dried(cmd,*args,**kwargs):
        #print, log and exit
        starttime_str=time.strftime("%y-%m-%d %H:%M:%S", time.localtime(time.time()))
        #convert cmd to string
        cmd=parse_cmd(cmd)
        command_name=cmd.split(' ')[0]
        log_message=cmd
        pu.print_blue("$ "+log_message)
        #log
        #create a dict and dump as json
        logDict={'cmd':log_message,
             'exitcode':"0",
             'runtime':"0",
             'starttime':str(starttime_str),
             'stdout':"dryrun",
             'stderr':"",
             'objectid':'',
             'commandname':command_name
            }
        if _logging:
            pyrpipe_logger.cmd_logger.debug(json.dumps(logDict))
        return True
    
    return dried

def parse_cmd(cmd):
    """This function converts a list to str.
    If a command is passed as list it is converted to str.
    pyrpipe v0.0.5 onwards the get_shell_output function uses shell=True
    """
    if isinstance(cmd,list):
        return " ".join(cmd)
    
    return cmd
        
@dryable
def get_shell_output(cmd,verbose=None):
    """Function to run a shell command and return returncode, stdout and stderr
    Currently (pyrpipe v 0.0.4) this function is called in 
    get_return_status(), get_program_version()
    
    pyrpipe v0.0.5 onwards the get_shell_output function allows shell=True
    
    
    Parameters
    ----------
    
    cdm: list or string
        command to run
    verbose: bool
        to print messages
    
    :return: (returncode, stdout and stderr)
    :rtype: tuple: (int,str,str)
    """
    if verbose==None: verbose=_verbose
    #not logging these commands
    cmd=parse_cmd(cmd)
    log_message=cmd
    command_name=cmd.split()[0]
    if _safe and command_name in ['rm']:
        pu.print_warning('SAFE MODE: Skipping command {}'.format(cmd))
        return True
    
    
    starttime_str=time.strftime("%y-%m-%d %H:%M:%S", time.localtime(time.time()))
    if verbose:
        pu.print_notification("Start:"+starttime_str)
        pu.print_blue("$ "+log_message)
    try:
        result = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
        stdout,stderr = result.communicate()
        if stdout:
            stdout=stdout.decode("utf-8")
        else:
            stdout=''
        if stderr:
            stderr=stderr.decode("utf-8")
        else:
            stderr=''
        return(result.returncode,stdout,stderr)
    except:
        return(-1,"Command failed to execute","Command failed to execute")

@dryable
def get_return_status(cmd):
    """
    run a shell command and get the return status
    
    Parameters
    ----------
    
    cmd: list
        shell command in list

    :return: True is returncode is 0
    :rtype: bool
    """
    #not logging these commands
    status=get_shell_output(cmd)
    if status[0]==0:
        return True
    return False

#prints stdout in real time. optimal for huge stdout and no stderr
@dryable
def execute_commandRealtime(cmd):
    """Execute shell command and print stdout in realtime.
    
    Example:
    for output in pe.execute_commandRealtime(['ping','-c','4','google.com']):
        print (output)
    """
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    
    for stdout_line in iter(popen.stdout.readline, ""):
            yield stdout_line 
    popen.stdout.close()
    return_code = popen.wait()
    
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

@dryable
def execute_command(cmd,verbose=None,logs=None,objectid=None,command_name=""):
    """Function to execute commands using popen. 
    All commands executed by this function can be logged and saved to pyrpipe logs.
    
    Parameters
    ----------
    
    cmd: list
        command to execute via popen in a list
    verbose: bool
        Whether to print stdout and stderr. Default: False. All stdout and stderr will be saved to logs regardless of this flag.
    logs: bool
        Log the execution 
    dryrun: bool
        If True, perform a dry run i.e. print commands to screen and log and exit
    objectid: string
        An id to be attached with the command. This is useful fo storing logs for SRA objects where object id is the SRR id.
    command_name: string
        Name of command to be save in log. If empty it is determined as the first element of the cmd list.

    :return: Return status.True is returncode is 0
    :rtype: bool
    """
    
    #if none then use default
    if verbose==None: verbose=_verbose
    if logs==None: logs=_logging
        
    
    #get current time
    time_start = time.time()
    starttime_str=time.strftime("%y-%m-%d %H:%M:%S", time.localtime(time.time()))
    
    if not command_name:
        command_name=cmd[0]
    log_message=" ".join(cmd)
    
    ###safe mode
    if _safe and command_name in ['rm']:
        pu.print_warning('SAFE MODE: Skipping command {}'.format(cmd))
        return True
    
    
    pu.print_notification("Start:"+starttime_str)
    pu.print_blue("$ "+log_message)
    
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
        
        endtime=time.time()
        endtime_str=time.strftime("%y-%m-%d %H:%M:%S", time.localtime(time.time()))
        timeDiff = round(endtime - time_start) #round to remove microsecond term
    
        if verbose:
            if stdout:
                pu.print_blue("STDOUT:\n"+stdout)
            if stderr:
                pu.print_boldred("STDERR:\n"+stderr)
        
        pu.print_notification("End:"+endtime_str)
        pu.print_green("Time taken:"+str(timedelta(seconds=timeDiff)))
            
                
        exitCode=result.returncode
        
        ##Add to logs        
        if logs:

            ##get the program used and log its path
            if command_name not in pyrpipe_logger.logged_programs:
                ##get which thisProgram
                #if subcommands are present use parent command
                parent_command=cmd[0]
                progDesc={'name':command_name,
                          'version':get_program_version(parent_command).strip(),
                          'path':get_program_path(parent_command).strip()
                          }
                pyrpipe_logger.env_logger.debug(json.dumps(progDesc))
                pyrpipe_logger.logged_programs.append(command_name)
            
            #create a dict and dump as json
            logDict={'cmd':log_message,
                 'exitcode':exitCode,
                 'runtime':str(timedelta(seconds=timeDiff)),
                 'starttime':str(starttime_str),
                 'stdout':stdout,
                 'stderr':stderr,
                 'objectid':str(objectid),
                 'commandname':command_name
                }
            pyrpipe_logger.cmd_logger.debug(json.dumps(logDict))
    
        if exitCode==0:
            return True
        else:
            #print the output
            pu.print_boldred("Following error occured executing above command (return code={}):".format(str(exitCode)))
            pu.print_boldred("STDOUT:\n"+stdout)
            pu.print_boldred("STDERR:\n"+stderr)
            return False
    #handle exceptions
    except OSError as e:
        pu.print_boldred("OSError exception occured.\n"+str(e))
        #log error
        timeDiff = round(time.time() - time_start)
        logDict={'cmd':log_message,
                 'exitcode':'-1',
                 'runtime':str(timedelta(seconds=timeDiff)),
                 'starttime':str(starttime_str),
                 'stdout':"",
                 'stderr':"OSError exception occured.\n"+str(e),
                 'objectid':objectid,
                 'commandname':command_name                 
                }
        pyrpipe_logger.cmd_logger.debug(json.dumps(logDict))
        return False
    except subprocess.CalledProcessError as e:
        pu.print_boldred("CalledProcessError exception occured.\n"+str(e))
        #log error
        timeDiff = round(time.time() - time_start)
        logDict={'cmd':log_message,
                 'exitcode':'-1',
                 'runtime':str(timedelta(seconds=timeDiff)),
                 'starttime':str(starttime_str),
                 'stdout':"",
                 'stderr':"CalledProcessError exception occured.\n"+str(e),
                 'objectid':objectid,
                 'commandname':command_name                 
                }
        pyrpipe_logger.cmd_logger.debug(json.dumps(logDict))
        return False
    except:
        pu.print_boldred("Fatal error occured during execution.\n"+str(sys.exc_info()[0]))
        #log error
        timeDiff = round(time.time() - time_start)
        logDict={'cmd':log_message,
                 'exitcode':'-1',
                 'runtime':str(timedelta(seconds=timeDiff)),
                 'starttime':str(starttime_str),
                 'stdout':"",
                 'stderr':str("Fatal error occured during execution.\n"+str(sys.exc_info()[0])),
                 'objectid':objectid,
                 'commandname':command_name
                }
        pyrpipe_logger.cmd_logger.debug(json.dumps(logDict))
        return False
    



def is_paired(sra_file):
    """Function to test wheather a .sra file is paired or single.
    
    Parameters
    ----------
    
    sra_file (string)  the path ro sra file
    
    
    :return: True is sra is paired
    :rtype: bool
    """
    if not pu.check_files_exist(sra_file):
        raise ValueError("Error checking layout. {0} doesn't exist".format(sra_file));
    
    try:
        fastqdCmd=["fastq-dump","-X","1","-Z","--split-spot", sra_file]
        if _dryrun:
            #if dry run 
            #pu.print_blue(' '.join(fastqdCmd))
            return True
        output = subprocess.check_output(fastqdCmd,stderr=subprocess.DEVNULL);
        numLines=output.decode("utf-8").count("\n")
        if(numLines == 4):
            return False;
        elif(numLines == 8):
            return True
        else:
            raise ValueError("Unexpected output from fast-dump");
    except subprocess.CalledProcessError as e:
        raise OSError("Error running fastq-dump: {}".format(str(e)));


    

def get_program_path(program):
    """
    Get path of installed program
    Returns the path as string    
    """
    return which(program)

def get_program_version(programName):
    """
    Get version of installed program
    return version as string
    """
    versionCommands=['--version','-version','--ver','-ver','-v','--v']
    for v in versionCommands:
        cmd=[programName,v]
        out=get_shell_output(cmd,verbose=False)
        if out[0]==0:
            return out[1]
    
    return ""

def check_dependencies(dependencies):
    """Check whether specified programs exist in the environment.
    This uses the which command to test whether a program is present.
    
    Parameters
    ----------
    
    dependencies: list
        list of programs to test

    :return: True is all dependencies are satified, False otherwise.
    :rtype: bool
    """
    for tool in dependencies:
        if not which(tool):
            pu.print_boldred("ERROR. {} command not found".format(tool))
            return False
    return True

  
    
#TODO: Re-implement following using native python libraries and move to utils
@skippable
def delete_file(file_path):
    """
    Delete a given file from disk
    Returns true if file is deleted or doesn't exist.
    shell=True is not added to make this function secure.

    Parameters
    ----------
    file_path : String
        Path to the file to be deleted
    Returns
    -------
    bool
        True if file deleted.

    """
    #if none
    if not file_path:
        return True
    
    #if pu.check_files_exist(file_path):
    rm_cmd=['rm',file_path]
    
    #rv= get_return_status(rm_cmd)
    rv=execute_command(rm_cmd)
    
    return rv
    
    #if file doesn't exist return true
    #return True

@skippable
def delete_files(*args):
    """Delete multiple files passed as argument.
    returns true is all files a re deleted
    """
    errorFlag=False
    for filePath in args:
        status=delete_file(filePath)
        if not status:
            errorFlag=True
    
    return not(errorFlag)

def move_file(source,destination,verbose=False):
    """perform mv command to move a file from sourc to destination
    Returns True if move is successful
    """
    if verbose:
        print("MOVING:"+source+"-->"+destination)
    mv_cmd=['mv',source,destination]
    if not get_return_status(mv_cmd):
        return False
    return True
