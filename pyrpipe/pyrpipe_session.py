#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 13:40:01 2019

@author: usingh
"""
import dill
import datetime as dt
import os

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
     
    
    #save workspace
    try:
        dill.dump_session(outFile)
        print("Session saved to: "+outFile)
        return True
    except Exception:
        raise OSError("Failed to save session")


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