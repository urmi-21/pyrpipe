#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 13:40:01 2019

@author: usingh
"""
import dill
import datetime as dt
import os
#importing pyrpipe_engine here causes issues and stalling of log
#from pyrpipe import pyrpipe_engine as pre
#import pyrpipe

def getTimestamp(shorten=False):
    """Return timestamp YYYYMMDDHHMISE
    shorten: return shorter version without spaces.
    """
    timestamp=str(dt.datetime.now()).split(".")[0].replace(" ","-")
    if shorten:
        timestamp=timestamp.replace("-","").replace(" ","").replace(":","")
    return timestamp

def save_session(filename,add_timestamp=True,out_dir=""):
    """Save current workspace using dill.
    Returns True is save is successful
    """
    #timestamp format YYYYMMDDHHMISE
    timestamp=getTimestamp(True)
    
    
    if not out_dir:        
        out_dir=os.getcwd()
    else:
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
    
    outFile=os.path.join(out_dir,filename)
    if add_timestamp:
        outFile=outFile+"_"+timestamp+".pyrpipe"
    else:
        outFile=outFile+".pyrpipe"
    
    """
    Do not pickle logger. This causes problems when restoring session with python < 3.7
    Delete all logger instances. 
    del pre.pyrpipeLoggerObject
    del pyrpipe.pyrpipe_engine.pyrpipeLoggerObject
    """    
    """
    creating a logger class fixed this issue
    """ 
    
    
    #save workspace
    try:
        dill.dump_session(outFile)
        print("Session saved to: "+outFile)
        return True
    except Exception:
        raise Exception("Failed to save session")
        return False


def restore_session(file):
    """Resore a session from file.
    
    :return: Returns True if session is restored
    :rtype: bool
    """
    if not os.path.isfile(file):
        print(file+" doesn't exist")
        return False
    #load the session
    dill.load_session(file)
    print("Session restored.")
    return True