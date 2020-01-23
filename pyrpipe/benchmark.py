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

class Benchmark:
    """Class to generate benchmark reports from pyrpipe logs.
    
    Parameters
    ----------
    
    log_file: string
        path to the log file
    env_log: string
        path to the ENV log file
    out_dir: string
        path to the output directory
        
        
    """
    def __init__(self,log_file,env_log,out_dir=""):
        
        
        if not pu.check_files_exist(log_file,env_log):
            raise Exception("Please check input for benchmark report. {} {}".format(log_file,env_log))
        if not out_dir:
            out_dir=os.getcwd()
        self.log_file=log_file
        self.env_log=env_log
        self.runtimes_by_prog={}
        self.runtimes_by_object={}
        #init
        pu.print_blue("parsing log...")
        self.parse_logs()
        pu.print_blue("done.")
        #out_dir
        self.benchmark_dir=os.path.join(out_dir,'benchmark_reports')
        if not pu.check_paths_exist(self.benchmark_dir):
            if not pu.mkdir(self.benchmark_dir):
                raise Exception("Error running benchmarks. Can not create output directory {}".format(self.benchmark_dir))
            
        
        
    def parse_runtime(self,timestring):
        """
        Parse runtime as string and return seconds.
        
        Returns: float
            runtime in sec
        """
        try:
            runtime= dt.datetime.strptime(timestring,"%H:%M:%S")
            delta_time = dt.timedelta(days=0,hours=runtime.hour, minutes=runtime.minute, seconds=runtime.second)
            return delta_time.seconds
        except ValueError:
            #try days format
            temp=timestring.split(",")
            days=int(temp[0].split(" ")[0].strip())
            rest=temp[1].strip()
            #hours=int(days)*24
            runtime= dt.datetime.strptime(rest,"%H:%M:%S")
            #one day less
            #lastruntime=lastruntime+dt.timedelta(days=days-1)
            delta_time = dt.timedelta(days=days,hours=runtime.hour, minutes=runtime.minute, seconds=runtime.second)
            return delta_time.seconds+(86400*days)
        
    def parse_logs(self):
        """Parse the input logs.
        For each command create a dict with runtimes as list and program name as key
        For each object id in input log create a dict.
        The followind dicts are created:
        
        runtimes_by_prog contains runtimes for each program. program is the key and the runtimes are in a list in order as they apprear in the log file.
        
        runtimes_by_object is nested a dict containing runtimes for each object by each program. e.g. {'ob1':{'prog1':[1,2,3],'prog2':[1,2,3]}, 'ob2':{'prog1':[12,22,13],'prog2':[1,2,3]} }
        """
        
        
        
        with open(self.log_file) as f:
            data=f.read().splitlines()
        num_commands=0

    
        for l in data:
            if not l.startswith("#"):
                thisDict=json.loads(l)
                num_commands+=1
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
                    objectid='SRR'+str(num_commands%50)
                
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
        """Returns a dataframe containing total execution time for each object in a pyrpipe log.
        An object is identified by the objectid e.g. SRR accession.
        """
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
        """Function to plot charts summarizing runtimes for each object in the pipeline.
        The charts are save to the out_dir path.
        """
        data=self.get_time_perobject()
        #remove rows with no object id
        data=data[data['id']!='NA']
        cols=list(data.columns)
        id_col=cols.index('id')
        #swap id with 0
        temp=cols[0]
        cols[0]='id'
        cols[id_col]=temp
        total_col=cols.index('total')
        #swap id with last
        temp=cols[-1]
        cols[-1]='total'
        cols[total_col]=temp
        
        data=data[cols]
        
        total_rows=data.shape[0]
        sns.set_context('poster')
        f, ax = plt.subplots(figsize = (total_rows/4,total_rows/2))
        
        sns.set_color_codes('bright')
        sns.barplot(x = 'total', y = 'id', data = data, label = 'Total', color = 'black', edgecolor = 'w')
        
        
        
        
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
        plotfile=os.path.join(self.benchmark_dir,'time_per_object.png')
        plt.savefig(plotfile,bbox_inches='tight')
        
        #write data to out_dir
        outfile=os.path.join(self.benchmark_dir,'time_per_object.csv')
        data.to_csv(outfile, index=False)
    
    
    def get_time_perprogram(self):
        """Returns a dataframe with program execution times.
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
    
    def get_programtime_boxdata(self):
        """Return dataframe to make box plot of program times
        """
        box_data=pd.DataFrame({ key:pd.Series(value) for key, value in self.runtimes_by_prog.items() })
        #return box_data
        
        
        ndf=pd.DataFrame(columns=['data','name'])
        for i in range(box_data.shape[1]):
            to_append=list(box_data.iloc[:,i])
            col_name=box_data.columns[i]
            col_list=[col_name for x in to_append]
            for j in range(len(to_append)):
                ndf=ndf.append({'data':to_append[j],'name':col_list[j]},ignore_index=True )
        
        return ndf
        
        
    def plot_time_perprogram(self):
        """Function to plot charts to summarize runtimes of each program
        """
        data=self.get_time_perprogram()
        sns.set_context('poster')
        f, ax = plt.subplots()
        
        current_palette = sns.color_palette("colorblind")
        
        #plot total time
        sns.barplot(x = 'total', y = 'program', data = data, label = 'Total', color = current_palette[4])
        #plot mean time
        sns.barplot(x = 'average', y = 'program', data = data, label = 'Avg.', color = current_palette[2])
        
        #save the barplot
        #add legend
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        #ax.legend(ncol = 2, loc = 'upper left')
        sns.despine(left = True, bottom = True)
        ax.set(xlabel='runtime (sec.)')
        #save plot
        plotfile=os.path.join(self.benchmark_dir,'time_per_program.png')
        plt.savefig(plotfile,bbox_inches='tight')
        
        
        #clear plot
        plt.clf()
        #make boxplots for runtime distribution
        sns.set_color_codes('bright')
        sns.set_context('poster')
        box_data=self.get_programtime_boxdata()
        #convert data to floats for boxplot
        box_data['data']=box_data['data'].astype(float)
        #print(box_data)
        numprog=len(self.runtimes_by_prog.keys())
        #sns.set(style="ticks")
        # Initialize the figure with a logarithmic x axis
        f, ax = plt.subplots(figsize=(20, numprog*2))
        ax.set_xscale("log")
        
        #sns.boxplot( data=box_data,orient="h")
        sns.boxplot(x="data", y="name", data=box_data)
        # Add in points to show each observation
        sns.swarmplot(x="data", y="name", data=box_data,size=5, color=".3", linewidth=0)
        
        ax.xaxis.grid(True)
        ax.set(ylabel="")
        ax.set(xlabel="runtime (sec.)")
        sns.despine(left=True)
        f.suptitle('Boxplot showing runtimes of each program', fontsize=26)
        
        
        #save to file
        plotfile=os.path.join(self.benchmark_dir,'program_boxplots.png')
        plt.savefig(plotfile,bbox_inches='tight')
        #plt.show()
        #clear plot
        
        plt.clf()
        
        
        #make pie charts
        sns.set_color_codes('bright')
        sns.set_context('poster')
        f, ax = plt.subplots(figsize=(10, numprog*2))
        current_palette = sns.color_palette("colorblind")
        plt.pie(data['total'], colors=current_palette, labels= data['program'],counterclock=False, shadow=False,autopct='%1.1f%%')
        f.suptitle('Summary of total runtimes', fontsize=26)
        
        #save plot
        plotfile=os.path.join(self.benchmark_dir,'program_summary.png')
        plt.savefig(plotfile,bbox_inches='tight')
        
        #write data to out_dir
        outfile=os.path.join(self.benchmark_dir,'time_per_program.csv')
        data.to_csv(outfile, index=False)
        #save boxplot data
        box_data=pd.DataFrame({ key:pd.Series(value) for key, value in self.runtimes_by_prog.items() })
        outfile=os.path.join(self.benchmark_dir,'program_box_data.csv')
        box_data.to_csv(outfile, index=False)
        
        
            
            


