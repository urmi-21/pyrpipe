#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 17:48:00 2019

@author: usingh
"""

from myutils import *

class Trimgalore:
    def __init__(self):
        self.programName="trim_galore"
        self.depList=[self.programName,'cutadapt']
        #check if hisat2 exists
        if not checkDep(self.depList):
            raise Exception("ERROR: "+ self.programName+" not found.")
            
    def runTrimGalorePaired(filepath,accession,proc):
        print ("Running trimgalore")
        #file names will be accession_1.fastq accession_2.fastq
        trimGaloreCmd=['trim_galore','-o',filepath,'--cores',str(proc),'--paired',filepath+"/"+accession+'_1.fastq',filepath+"/"+accession+'_2.fastq']
        print("Executing: "+" ".join(trimGaloreCmd))
        try:
            	for output in executeCommand(trimGaloreCmd):
                    nt (output)
        except subprocess.CalledProcessError as e:
            	print ("Error in command")
            return False
    print("Exiting...")
    return True
            
            

class BBmap:
    def __init__(self):
        self.programName="bbduk.sh"
        self.depList=[self.programName]
        #check if hisat2 exists
        if not checkDep(self.depList):
            raise Exception("ERROR: "+ self.programName+" not found.")
            
            
    def runBBDUK(filepath,accession,pathToAdapters,proc="auto",ktrim='r',k=23,mink=11,hdist=1,qtrim='rl',trimq=10):
        print ("Running bbduk")
        #file names will be accession_1.fastq accession_2.fastq
        bbdukCmd=['bbduk.sh','-Xmx1g','in1='+filepath+"/"+accession+'_1.fastq','in2='+filepath+"/"+accession+'_2.fastq','out1='+filepath+"/"+accession+'_1_val_1.fastq','out2='+filepath+"/"+accession+'_2_val_2.fastq','ref='+pathToAdapters,'ktrim='+ktrim,'k='+str(k),'mink='+str(mink),'hdist='+str(hdist),'qtrim='+qtrim,'trimq='+str(trimq),'threads='+str(proc)]
        print("Executing: "+" ".join(bbdukCmd))
        try:
            	for output in executeCommand(bbdukCmd):
                    print (output)
        except subprocess.CalledProcessError as e:
            	print ("Error in command")
            return False
        print("Exiting...")
        return True