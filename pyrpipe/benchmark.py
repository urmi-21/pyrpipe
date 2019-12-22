#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 16:58:13 2019

@author: usingh
"""

from pyrpipe import pyrpipe_utils as pu
import datetime as dt
import json
import seaborn as sns
import pandas as pd
import copy 
import matplotlib.pyplot as plt


class benchmark:
    def __init__(self,logFile,envLog):
        
        if not pu.checkFilesExists(logFile,envLog):
            raise Exception("Please check input for benchmark report. {} {}".format(logFile,envLog))
        self.logFile=logFile
        self.envLpg=envLog
        self.parseLogs()
        
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
        
    def parseLogs(self):
        """Parse the input logs.
        For each command create a dict with runtimes as list and program name as key
        For each object id in input log create a dict.
        """
        
        self.runtimes_by_prog={}
        self.runtimes_by_object={}
        
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
                             
                #store runtimes by programname
                programname=thisDict['commandname']
                runtime=self.parse_runtime(thisDict['runtime'])
                #add to dict
                if programname in self.runtimes_by_prog:
                    self.runtimes_by_prog[programname].append(runtime)
                else:
                    self.runtimes_by_prog[programname]=[runtime]
                    
                #store runtimes by object id
                objectid=thisDict['objectid']
                
                if objectid in self.runtimes_by_object:
                    #if the object is used with same program extend the list
                    if programname in self.runtimes_by_object[objectid]:
                        self.runtimes_by_object[objectid][programname].append(runtime)
                    else:
                        self.runtimes_by_object[objectid][programname]=[runtime]
                else:
                    self.runtimes_by_object[objectid]={programname:[runtime]}
                
        #print(self.runtimes_by_prog)
        #print(self.runtimes_by_object)
                                   
    def get_time_perobject(self,func="sum"):
        result=pd.DataFrame()
        for k in self.runtimes_by_object:
            
            v=copy.deepcopy(self.runtimes_by_object[k]) #don't change the self.runtimes_by_object
            
            totaltime=0
            for prog in v:
                runtimes=v[prog]                
                if func == "mean":
                    valtoput=sum(runtimes)/len(runtimes)
                    totaltime+=valtoput
                    v[prog]=[valtoput]
                else:                    
                    valtoput=sum(runtimes)
                    totaltime+=valtoput
                    v[prog]=[valtoput]
                    
            #add current id
            v['id']=k
            v['total']=totaltime
            
            result=result.append(pd.DataFrame.from_dict(v),sort=False)
            
        return result
        
    def plot_time_perobject(self):
        data=self.get_time_perobject()
        sns.set_context('paper')
        f, ax = plt.subplots(figsize = (6,15))
        
        sns.set_color_codes('bright')
        sns.barplot(x = 'total', y = 'id', data = data, label = 'Total', color = 'black', edgecolor = 'w')
        
        
        #sns.set_color_codes('muted')
        
        #for each col
        current_palette = sns.color_palette("colorblind")
        i=0
        for col in data.columns:
            if col in ['id','total']:
                continue
            sns.barplot(x = col, y = 'id', data = data,label = col, color = current_palette[i], edgecolor = 'w',)
            i+=1

        ax.legend(ncol = 2, loc = 'lower right')
        sns.despine(left = True, bottom = True)
        plt.show()
            
            
if __name__ == "__main__":
    print("testing")
    l="/home/usingh/work/urmi/hoap/pyrpipe/tests/pyrpipe_logs/2019-12-21-16_38_08_pyrpipe.log"
    e="/home/usingh/work/urmi/hoap/pyrpipe/tests/pyrpipe_logs/2019-12-21-16_38_08_pyrpipeENV.log"
    ob=benchmark(l,e)
    d=(ob.get_time_perobject())
    ob.plot_time_perobject()







