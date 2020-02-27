#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 16:26:02 2019

@author: usingh

Methods and classes related to execution and logging
"""

import os
import subprocess
import time
from datetime import datetime 
from datetime import timedelta
import logging
import sys
import platform
from multiprocessing import cpu_count
from pyrpipe import pyrpipe_utils as pu
import json

class LogFormatter():
    """
    A formatter for logs
    """
    def __init__(self):
        self.start_time = time.time()
    
    def format(self, record):
      
        return "{}".format(record.getMessage())       
        #time_now=str(datetime.now())
        #elapsed_seconds = record.created - self.start_time
        #hours, remainder = divmod(elapsed_seconds, 3600)
        #minutes, seconds = divmod(remainder, 60)
        
        #message=record.getMessage()

        """OLD
        if record.name == "env":
            return "{}".format(record.getMessage())
        
        if message.startswith("#"):
            return "{}\t{}".format(record.getMessage(),time_now)
        """
        
        """
        if record.name == "cmd":
            return "{}\nStart:{}\tDuration: {:02}:{:02}:{:02}".format(record.getMessage(),time_now, int(hours), int(minutes), int(seconds) )
        elif record.name == "out":
            return "STDOUT:\n{}".format(record.getMessage())
        elif record.name == "err":      
            return "STDERR:\n{}".format(record.getMessage())
        """
        

class PyrpipeLogger():
    """
    Class to manage pyrpipe logs
    
    Attributes
    -----------
    
    env_logger: logger to log the current environment
    cmd_logger: logger to log the execution status, stdout, stderr and runtimes for each command run using execute_command()
    """
    def __init__(self):
        self.__name__="pyrpipeLogger"
        #loggers
        timestamp=str(datetime.now()).split(".")[0].replace(" ","-").replace(":","_")
        self.logger_basename=timestamp+"_pyrpipe"
        self.logs_dir=os.path.join(os.getcwd(),"pyrpipe_logs")
        if not os.path.isdir(self.logs_dir):
            os.mkdir(self.logs_dir)
        """
        self.cmd_loggerPath=os.path.join(self.logs_dir,self.logger_basename+"CMD.log")
        self.stdoutLoggerPath=os.path.join(self.logs_dir,self.logger_basename+"OUT.log")
        self.stderrLoggerPath=os.path.join(self.logs_dir,self.logger_basename+"ERR.log")
        """
        self.log_path=os.path.join(self.logs_dir,self.logger_basename+".log")
        self.envlog_path=os.path.join(self.logs_dir,self.logger_basename+"ENV.log")
        
        """
        self.cmd_logger=self.create_logger("cmd",self.cmd_loggerPath,LogFormatter(),logging.DEBUG)
        self.stdoutLogger=self.create_logger("out",self.stdoutLoggerPath,LogFormatter(),logging.DEBUG)
        self.stderrLogger=self.create_logger("err",self.stderrLoggerPath,LogFormatter(),logging.DEBUG)
        """
        formatter=LogFormatter()
        self.env_logger=self.create_logger("env",self.envlog_path,formatter,logging.DEBUG)
        self.cmd_logger=self.create_logger("cmd",self.log_path,formatter,logging.DEBUG)
        
        #self.stdoutLogger=self.create_logger("out",self.log_path,formatter,logging.DEBUG)
        #self.stderrLogger=self.create_logger("err",self.log_path,formatter,logging.DEBUG)
        
        #self.env_logger=self.create_logger("env",self.envlog_path,logging.Formatter("%(message)s"),logging.DEBUG)
        
        self.init_envlog()
        self.init_cmdlog()
 
        
    
    def create_logger(self,name,logfile,formatter,level=logging.DEBUG):
        """Creates a logger
        
        Parameters
        ----------
        
        name: str
            name of logger
        logfile: str
            file name to save logs
        formatter: formatter object
            formatter for log
        
        Returns: logger
            A logger object
        """
        #Get different loggers
        handler = logging.FileHandler(logfile)        
        handler.setFormatter(formatter)
        
        logger = logging.getLogger(name)
        logger.setLevel(level)
        logger.addHandler(handler)
        return logger
    
    def init_cmdlog(self):
        """init the cmdlog
        """
        self.cmd_logger.debug("#START LOG")
    

    def init_envlog(self):
        """init the envlog
        """
        self.env_logger.debug("#START LOG")
        
        #get current time
        #sesstime='Session information collected on {}'.format(datetime.now().strftime('%Y-%m-%d %H:%M'))
        #get os
        osInfo=platform.platform()
        #get python version
        pyver='Python ' + sys.version.replace('\n', '')
        #get cpu
        cpu=str(cpu_count())+' logical CPU cores'
        
        envDesc={'now':str(datetime.now().strftime("%y-%m-%d %H:%M:%S")),
                 'python':pyver,
                 'os':osInfo,
                 'cpu':cpu,
                 'syspath':str(sys.path),
                 'sysmodules':str(list(sys.modules.keys()))
                 }
        
        self.env_logger.debug(json.dumps(envDesc))
        """Old
        self.env_logger.debug(sesstime)
        self.env_logger.debug(pyver)
        self.env_logger.debug(osInfo)
        self.env_logger.debug(cpu)
        self.env_logger.debug("#SYS PATH")
        self.env_logger.debug("sys.path:"+str(sys.path))
        self.env_logger.debug("#SYS MODULES")
        self.env_logger.debug("sys.modules:"+str(sys.modules.keys()))
        """
        self.env_logger.debug("#PROGRAMS")
        #a list of logged programs
        self.logged_programs=[]
        

###create logger
pyrpipeLoggerObject=PyrpipeLogger()
pu.print_yellow("Logs will be saved to {}.log".format(pyrpipeLoggerObject.logger_basename))
    
"""
All functions that interact with shell are defined here. 
"""
    
##############Functions###########################

   
    
def getShellOutput(cmd,verbose=False):
    """Function to run a shell command and return returncode, stdout and stderr
    
    Parameters
    ----------
    
    cdm: list
        command to run
    verbose: bool
        to print messages
    
    :return: (returncode, stdout and stderr)
    :rtype: tuple
    """
    #not logging these commands
    log_message=" ".join(cmd)
    if verbose:
        pu.print_blue("$ "+log_message)
    try:
        result = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,)
        stdout,stderr = result.communicate()
        return(result.returncode,stdout,stderr)
    except:
        return(-1,"","Command failed to execute")

def getReturnStatus(cmd):
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
    status=getShellOutput(cmd)
    if status[0]==0:
        return True
    return False

#prints stdout in real time. optimal for huge stdout and no stderr
def execute_commandRealtime(cmd):
    """Execute shell command and print stdout in realtime.
    """
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    
    for stdout_line in iter(popen.stdout.readline, ""):
            yield stdout_line 
    popen.stdout.close()
    return_code = popen.wait()
    
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


def execute_command(cmd,verbose=False,quiet=False,logs=True,dryrun=False,objectid="NA",command_name=""):
    """Function to execute commands using popen. 
    All commands executed by this function can be logged and saved to pyrpipe logs.
    
    Parameters
    ----------
    
    cmd: list
        command to execute via popen in a list
    verbose: bool
        Whether to print stdout and stderr. Default: False. All stdout and stderr will be saved to logs regardless of this flag.
    quiet: bool
        Absolutely no output on screen
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
    if not command_name:
        command_name=cmd[0]
    log_message=" ".join(cmd)
    
    #dryrun: print and exit
    if dryrun:
        pu.print_blue("$ "+log_message)
        #log
        #create a dict and dump as json
        logDict={'cmd':log_message,
                 'exitcode':"0",
                 'runtime':"0",
                 'starttime':"0",
                 'stdout':"dryrun",
                 'stderr':"",
                 'objectid':objectid,
                 'commandname':command_name
                }
        pyrpipeLoggerObject.cmd_logger.debug(json.dumps(logDict))
        return True
    
    if not quiet:
        pu.print_blue("$ "+log_message)
    time_start = time.time()
    starttime_str=time.strftime("%y-%m-%d %H:%M:%S", time.localtime(time.time()))
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
        
        timeDiff = round(time.time() - time_start) #round to remove microsecond term
    
        if verbose:
            if stdout:
                pu.print_blue("STDOUT:\n"+stdout)
            if stderr:
                pu.print_boldred("STDERR:\n"+stderr)
        if not quiet:
            pu.print_green("Time taken:"+str(timedelta(seconds=timeDiff)))
            
                
        exitCode=result.returncode
        
        ##Add to logs        
        if logs:

            ##get the program used and log its path
            if command_name not in pyrpipeLoggerObject.logged_programs:
                ##get which thisProgram
                #if subcommands are present use parent command
                parent_command=cmd[0]
                progDesc={'name':command_name,
                          'version':getProgramVersion(parent_command).strip(),
                          'path':getProgramPath(parent_command).strip()
                          }
                pyrpipeLoggerObject.env_logger.debug(json.dumps(progDesc))
                pyrpipeLoggerObject.logged_programs.append(command_name)
            
            #create a dict and dump as json
            logDict={'cmd':log_message,
                 'exitcode':exitCode,
                 'runtime':str(timedelta(seconds=timeDiff)),
                 'starttime':str(starttime_str),
                 'stdout':stdout,
                 'stderr':stderr,
                 'objectid':objectid,
                 'commandname':command_name
                }
            pyrpipeLoggerObject.cmd_logger.debug(json.dumps(logDict))
    
        if exitCode==0:
            return True
        else:
            #print the output
            print("Following error occured executing above command (return code={}):".format(str(exitCode)))
            print("STDOUT:\n"+stdout)
            print("STDERR:\n"+stderr)
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
        pyrpipeLoggerObject.cmd_logger.debug(json.dumps(logDict))
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
        pyrpipeLoggerObject.cmd_logger.debug(json.dumps(logDict))
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
        pyrpipeLoggerObject.cmd_logger.debug(json.dumps(logDict))
        return False
    



#modified from https://www.biostars.org/p/139422/
def is_paired(sra_file):
    """Function to test wheather a .sra file is paired or single.
    
    Parameters
    ----------
    
    sra_file (string)  the path ro sra file
    
    
    :return: True is sra is paired
    :rtype: bool
    """
    if not pu.check_files_exist(sra_file):
        raise Exception("Error checking layout. {0} doesn't exist".format(sra_file));
    
    try:
        fastqdCmd=["fastq-dump","-X","1","-Z","--split-spot", sra_file]
        output = subprocess.check_output(fastqdCmd,stderr=subprocess.DEVNULL);
        numLines=output.decode("utf-8").count("\n")
        if(numLines == 4):
            return False;
        elif(numLines == 8):
            return True
        else:
            raise Exception("Unexpected output from fast-dump");
    except subprocess.CalledProcessError as e:
        raise Exception("Error running fastq-dump: {}".format(str(e)));


    

def getProgramPath(programName):
    """
    Get path of installed program
    Returns the path as string    
    """
    whichCmd=['which',programName]
    out = subprocess.check_output(whichCmd,universal_newlines=True)
    return out

def getProgramVersion(programName):
    """
    Get version of installed program
    return version as string
    """
    versionCommands=['--version','-version','--ver','-ver','-v','--v']
    for v in versionCommands:
        cmd=[programName,v]
        out=getShellOutput(cmd)
        if out[0]==0:
            return out[1].decode("utf-8")
    
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
    errorFlag=False
    for s in dependencies:
        #print_blue("Checking "+s+"...")
        thisCmd=['which',s]
        if(getReturnStatus(thisCmd)):
            #print_green ("Found "+s)
            pass
        else:
            pu.print_boldred ("Can not find "+s)
            errorFlag=True
    if errorFlag:
        return False
    return True


#TODO: Re-implement following using native python libraries and move to utils
def deleteFileFromDisk(filePath):
    """Delete a given file from disk
    Returns true if file is deleted or doesn't exist
    """
    if pu.check_files_exist(filePath):
        rm_Cmd=['rm',filePath]
        rv= getReturnStatus(rm_Cmd)
        return rv
    #if file doesn't exist return true
    return True

def deleteMultipleFilesFromDisk(*args):
    """Delete multiple files passed as argument.
    returns true is all files a re deleted
    """
    errorFlag=False
    for filePath in args:
        status=deleteFileFromDisk(filePath)
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
    if not getReturnStatus(mv_cmd):
        return False
    return True


#TODO: use os.scandir
def find_files(search_path,search_pattern,recursive=False,verbose=False):
    """Function to find files using find command and return as list
    Use global paths for safety
    
    Parameters
    ----------
    
    search_path: str
        path to search under
    search_pattern: str
        pattern to search e.g. "*.bam"
    recursive: bool
        search all subdirs if true

    :return: list containing the found paths
    :rtype: list
    """
    
    if recursive:
        find_cmd=['find', search_path,'-type', 'f','-name',search_pattern]   
    else:
        find_cmd=['find', search_path,'-maxdepth', '1','-type', 'f','-name',search_pattern] 
    if verbose:
        print("$ "+" ".join(find_cmd))
    st=getShellOutput(find_cmd)
    output=[]
    if st[0]==0:
        output=st[1].decode("utf-8").split("\n")
    #remove empty strings
    if '' in output:
        output.remove('')
    return output
    



