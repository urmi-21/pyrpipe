#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 16:58:13 2019

@author: usingh
"""

from pyrpipe import pyrpipe_utils as pu
import datetime as dt
class benchmark:
    def __init__(self,logFile,envLog):
        
        if not pu.checkFilesExists(logFile,envLog):
            raise Exception("Please check input for benchmark report. {} {}".format(logFile,envLog))
        self.logFile=logFile
        self.envLpg=envLog
        
        
        
    def parseLogs(self):
        """Parse the input logs.
        For each command create a dict with runtimes as list and program name as key
        For each object id in input log create a dict.
        """
        
        self.times_by_prog={}
        
        with open(self.logFile) as f:
        data=f.read().splitlines()
        numCommands=0
        failedCommands=0
        passedCommands=0
        numPrograms=0
    
        for l in data:
            if not l.startswith("#"):
            thisDict=json.loads(l)
            numCommands+=1
            #ignore failed commands
            if int(thisDict['exitcode'])!=0:
                continue
                                   
            #program name
            programname=thisDict['commandname']
            runtime=self.parse_runtime(thisDict['runtime'])
            #add to dict
            if programname in self.times_by_prog:
                self.times_by_prog[programname].append(runtime)
            else:
                self.times_by_prog[programname]=[runtime]
                                   
            
        def parse_runtime(self,timestring):
            """
            Parse runtime as string and return seconds.
            """
            try:
                runtime= dt.datetime.strptime(timestring,"%H:%M:%S")
                deltaTime = dt.timedelta(days=0,hours=runtime.hour, minutes=runtime.minute, seconds=runtime.second)
                return deltaTime.seconds
            except ValueError:
                #try days format
                temp=timestring.split(",")
                days=int(temp[0].split(" ")[0].strip())
                rest=temp[1].strip()
                #hours=int(days)*24
                runtime= dt.datetime.strptime(rest,"%H:%M:%S")
                #one day less
                #lastruntime=lastruntime+dt.timedelta(days=days-1)
                deltaTime = dt.timedelta(days=days,hours=runtime.hour, minutes=runtime.minute, seconds=runtime.second)
                return deltaTime.seconds+(86400*days)
            
            
            
            
if __name__ == "__main__":
    print("testing")
            
            
            







