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

import os

class benchmark:
    def __init__(self,logFile,envLog,outDir=os.getcwd()):
        """Describe...
        
        
        
        """
        
        if not pu.checkFilesExists(logFile,envLog):
            raise Exception("Please check input for benchmark report. {} {}".format(logFile,envLog))
        self.logFile=logFile
        self.envLpg=envLog
        self.runtimes_by_prog={}
        self.runtimes_by_object={}
        #init
        pu.printBlue("parsing log...")
        self.parseLogs()
        pu.printBlue("done.")
        #outdir
        self.benchmarksDir=os.path.join(outDir,'benchmark_reports')
        if not pu.checkPathsExists(self.benchmarksDir):
            if not pu.mkdir(self.benchmarksDir):
                raise Exception("Error running benchmarks. Can not create output directory {}".format(self.benchmarksDir))
            
        
        
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
        The followind dicts are created:
        
        runtimes_by_prog contains runtimes for each program. program is the key and the runtimes are in a list in order as they apprear in the log file.
        
        runtimes_by_object is nested a dict containing runtimes for each object by each program. e.g. {'ob1':{'prog1':[1,2,3],'prog2':[1,2,3]}, 'ob2':{'prog1':[12,22,13],'prog2':[1,2,3]} }
        """
        
        
        
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
                try:
                    programname=thisDict['commandname']
                except KeyError:
                    #for older logs
                    programname=thisDict['cmd'].split(" ")[0]
                runtime=self.parse_runtime(thisDict['runtime'])
                #add to dict
                if programname in self.runtimes_by_prog:
                    self.runtimes_by_prog[programname].append(runtime)
                else:
                    self.runtimes_by_prog[programname]=[runtime]
                    
                #store runtimes by object id
                try:
                    objectid=thisDict['objectid']
                except KeyError:
                    objectid='SRR'+str(numCommands%50)
                
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
        #remove rows with no object id
        data=data[data['id']!='NA']
        total_rows=data.shape[0]
        sns.set_context('poster')
        f, ax = plt.subplots(figsize = (total_rows/4,total_rows/2))
        
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

        #add legend
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        #ax.legend(ncol = 2, loc = 'upper right')
        sns.despine(left = True, bottom = True)
        ax.set(xlabel='runtime (sec.)')
        #plt.show()
        
        #save plot
        plotfile=os.path.join(self.benchmarksDir,'time_per_object.png')
        plt.savefig(plotfile,bbox_inches='tight')
        
        #write data to outdir
        outfile=os.path.join(self.benchmarksDir,'time_per_object.csv')
        data.to_csv(outfile)
    
    
    def get_time_perprogram(self):
        """Returns a dataframe with program execution times
        """
        result=pd.DataFrame()
        for k in self.runtimes_by_prog:
            v=self.runtimes_by_prog[k]
            total=sum(v)
            mean=sum(v)/len(v)
            row={'program':k,'total':[total],'average':[mean]}
            #add row to dataframe
            result=result.append(pd.DataFrame.from_dict(row),sort=False)
        return result
            
        
        
    def plot_time_perprogram(self):
        data=self.get_time_perprogram()
        sns.set_context('poster')
        f, ax = plt.subplots()
        
        sns.set_color_codes('bright')
        
        #plot total time
        sns.barplot(x = 'total', y = 'program', data = data, label = 'Total', color = 'black', edgecolor = 'black')
        #plot mean time
        sns.barplot(x = 'average', y = 'program', data = data, label = 'Avg.', color = 'red', edgecolor = 'red')
        
        #save the barplot
        #add legend
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        #ax.legend(ncol = 2, loc = 'upper left')
        sns.despine(left = True, bottom = True)
        ax.set(xlabel='runtime (sec.)')
        #save plot
        plotfile=os.path.join(self.benchmarksDir,'time_per_program.png')
        plt.savefig(plotfile,bbox_inches='tight')
        
        #clear plot
        plt.clf()
        
        
        #make pie charts
        current_palette = sns.color_palette("colorblind")
        plt.pie(data['total'], colors=current_palette, labels= data['program'],counterclock=False, shadow=True)
        #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        
        #save plot
        plotfile=os.path.join(self.benchmarksDir,'program_summary.png')
        plt.savefig(plotfile,bbox_inches='tight')
        
        #write data to outdir
        outfile=os.path.join(self.benchmarksDir,'time_per_program.csv')
        data.to_csv(outfile)
        
        
            
            
if __name__ == "__main__":
    print("testing")
    l="/home/usingh/work/urmi/hoap/test/bmtest/2019-12-19-13_17_11_pyrpipe.log"
    e="/home/usingh/work/urmi/hoap/test/bmtest/2019-12-19-13_17_11_pyrpipeENV.log"
    ob=benchmark(l,e,outDir="/home/usingh/work/urmi/hoap/test/bmtest")
    d=(ob.get_time_perobject())
    ob.plot_time_perobject()

    d2=ob.get_time_perprogram()
    ob.plot_time_perprogram()





