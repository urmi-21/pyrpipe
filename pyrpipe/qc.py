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
        self.trim_galoreArgsList=['-h','-v','-q','--phred33   ','--phred64','--fastqc ','--fastqc_args','-a','-a2','--illumina','--nextera','--small_rna','--consider_already_trimmed         ','--max_length','--stringency','-e','--gzip','--dont_gzip','--length','--max_n ','--trim-n','-o','--no_report_file','--suppress_warn','--clip_R1','--clip_R2','--three_prime_clip_R1','--three_prime_clip_R2','--2colour','--path_to_cutadapt','--basename','-j','--hardtrim5','--hardtrim3','--clock','--polyA','--rrbs','--non_directional','--keep','--paired','-t','--retain_unpaired','-r1','-r2']
        #check if hisat2 exists
        if not checkDep(self.depList):
            raise Exception("ERROR: "+ self.programName+" not found.")
            
            
    def runTrimGalore(self,sraOb):
        #get layout
        if sraOb.layout=='PAIRED':
            fq1=sraOb.localfastq1Path
            fq2=sraOb.localfastq2Path
            return self.runTrimGalorePaired(fq1,fq2,8)
        else:
            return self.runTrimGaloreSingle(sraOb.localfastqPath,8)
            
        
            
    def runTrimGaloreSingle(self,fastqFilePath,proc):
        
        #default output dir
        outDir=os.path.split(fastqFilePath)[0]
        
        print("Running trim_galore unpaired")
        trimGaloreCmd=['trim_galore','-o',outDir,'--cores',str(proc),fastqFilePath]
        
            
    def runTrimGalorePaired(self,fastqFile1Path,fastqFile2Path,proc):
        print ("Running trim_galore paired")
        #default output dir
        outDir=os.path.split(fastqFile1Path)[0]
        
        trimGaloreCmd=['trim_galore','-o',outDir,'--cores',str(proc),'--paired',fastqFile1Path,fastqFile1Path]
        print("Executing: "+" ".join(trimGaloreCmd))
        try:
            for output in executeCommand(trimGaloreCmd):
                print (output)
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
            
            
    def runBBDUK(self,filepath,accession,pathToAdapters,proc="auto",ktrim='r',k=23,mink=11,hdist=1,qtrim='rl',trimq=10):
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