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
    
    timestamp=str(dt.datetime.now()).split(".")[0].replace(" ","-")
    if shorten:
        timestamp=timestamp.replace("-","").replace(" ","").replace(":","")
    return timestamp

def save_session(filename="session",timestamp=True,out_dir=""):
    """Save current workspace using dill.
    """
    #timestamp format YYYYMMDDHHMISE
    timestamp=getTimestamp(True)
    
    
    if not out_dir:        
        out_dir=os.getcwd()
    else:
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
    
    outFile=os.path.join(out_dir,filename)
    if timestamp:
        outFile=outFile+"_"+timestamp+".pyrpipe"
    
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
    dill.dump_session(outFile)
    print("Session saved to: "+outFile)
    
    return True


def restore_session(file):
    if not os.path.isfile(file):
        print(file+" doesn't exist")
        return False
    #load the session
    dill.load_session(file)
    print("Session restored.")
    return True