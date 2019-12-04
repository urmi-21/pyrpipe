#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 17:48:00 2019

@author: usingh
"""

from pyrpipe.myutils import *

class RNASeqQC:
    def __init__(self):
        self.category="RNASeqQC"
    def performQC(self):
        pass

class Trimgalore(RNASeqQC):
    def __init__(self,**kwargs):
        """
        Parameters
        ----------
        kwargs:
            trim_galore arguments. could override later too.
        """
        
        #run super to inherit parent class properties
        super().__init__() 
        self.programName="trim_galore"
        self.depList=[self.programName,'cutadapt']
        self.validArgsList=['-h','-v','-q','--phred33','--phred64','--fastqc','--fastqc_args','-a','-a2',
                            '--illumina','--nextera','--small_rna','--consider_already_trimmed',
                            '--max_length','--stringency','-e','--gzip','--dont_gzip','--length',
                            '--max_n','--trim-n','-o','--no_report_file','--suppress_warn',
                            '--clip_R1','--clip_R2','--three_prime_clip_R1','--three_prime_clip_R2',
                            '--2colour','--path_to_cutadapt','--basename','-j','--hardtrim5','--hardtrim3',
                            '--clock','--polyA','--rrbs','--non_directional','--keep','--paired','-t',
                            '--retain_unpaired','-r1','-r2']
        #check if hisat2 exists
        if not checkDep(self.depList):
            raise Exception("ERROR: "+ self.programName+" not found.")
            
        #initialize the passed arguments
        self.passedArgumentDict=kwargs
        
        
            
    def performQC(self,sraOb,outFileSuffix="_trimgalore",**kwargs):
        """Function to perform qc using trimgalore.
        The function performQC() is consistent for all QC classess.
        
        Parameters
        ----------
        arg1: SRA
            An SRA object whose fastq files will be used
        
        Returns
        -------
        tuple
            Returns the path of fastq files after QC. tuple has one item for single end files and 2 for paired.
            
        """
        #create path to output files
        outDir=sraOb.location
        
        #create new options based on parametrs
        newOpts={}
        #get layout
        if sraOb.layout=='PAIRED':
            fq1=sraOb.localfastq1Path
            fq2=sraOb.localfastq2Path
            outFile1=os.path.join(outDir,getFileBaseName(fq1)+outFileSuffix+".fastq")
            outFile2=os.path.join(outDir,getFileBaseName(fq2)+outFileSuffix+".fastq")
            newOpts={"--paired":"","--":(fq1,fq2),"-o":outDir}
            mergedOpts={**kwargs,**newOpts}
            #run trimgalore
            self.runTrimGalore(**mergedOpts)
            """
            running trim galore will create two files named <input>_val_1.fq and <input>_val_2.fq
            move these files to the specified out files
            """
            oldFile1=os.path.join(outDir,getFileBaseName(fq1)+"_val_1.fq")
            oldFile2=os.path.join(outDir,getFileBaseName(fq2)+"_val_2.fq")
            
            mv1=moveFile(oldFile1,outFile1)
            mv2=moveFile(oldFile2,outFile2)
            
            if not checkFilesExists(outFile1,outFile2):
                print("Trimgalore failed")
                return ("",)
            return outFile1,outFile2
            
        else:
            fq=sraOb.localfastqPath
            outFile=os.path.join(outDir, getFileBaseName(fq)+outFileSuffix+".fastq")
            #giving input arguments as a tuple "--":(fq,)
            newOpts={"--":(fq,),"-o":outDir}
            #run trimgalore
            mergedOpts={**kwargs,**newOpts}
            print("RUN Single TG:"+str(mergedOpts))
            self.runTrimGalore(**mergedOpts)
            """
            running trim galore will create one file named <input>_trimmed.fq
            move these files to the specified out files
            """
            oldFile=os.path.join(outDir,getFileBaseName(fq)+"_trimmed.fq")
            
            mv=moveFile(oldFile,outFile)
            
            if not checkFilesExists(outFile):
                print("Trimgalore failed")
                return ("",)
            return (outFile,)
        
        
            
    def runTrimGalore(self,**kwargs):
        """Wrapper for running trimgalore
        
        Parameters
        ----------
        arg1: dict
            Options to pass to trimgalore (will override existing parameters)
        """
        
        #override existing arguments
        mergedArgsDict={**self.passedArgumentDict,**kwargs}
        
        #create command to run
        trimGalore_Cmd=['trim_galore']
        trimGalore_Cmd.extend(parseUnixStyleArgs(self.validArgsList,mergedArgsDict))
        print("Executing:"+" ".join(trimGalore_Cmd))
        
        #start ececution
        log=""
        try:
            for output in executeCommand(trimGalore_Cmd):
                #print (output)    
                log=log+str(output)
            #save to a log file
        except subprocess.CalledProcessError as e:
            print ("Error in command...\n"+str(e))
            #save error to error.log file
            return False        
        #return status
        return True
        

            
            

class BBmap(RNASeqQC):
    def __init__(self):
        """
        Parameters
        ----------
        kwargs:
            bbduk.sh arguments. could override later too.
        """
        #run super to inherit parent class properties
        super().__init__() 
        self.programName="bbduk.sh"
        self.depList=[self.programName]
        #note that bbduk.sh argument style is different that other linux commands
        self.validArgsList=[]
        #check if hisat2 exists
        if not checkDep(self.depList):
            raise Exception("ERROR: "+ self.programName+" not found.")
            
    def run(self,sraOb):
        """Execeute the QC method 
        """
        if sraOb.layout=='PAIRED':
            fq1=sraOb.localfastq1Path
            fq2=sraOb.localfastq2Path
            return self.runBBdukPaired(fq1,fq2,"/home/usingh/lib_urmi/softwares/bbmap/resources/adapters2.fa")
        else:
            return self.runBBdukSingle(sraOb.localfastqPath,"/home/usingh/lib_urmi/softwares/bbmap/resources/adapters2.fa")
    
    #default ktrim='r',k=23,mink=11,hdist=1,qtrim='rl',trimq=10
    def runBBdukSingle(self,fastqFilePath,pathToAdapters="",proc="auto",ktrim='r',k=13,mink=5,hdist=1,qtrim='rl',trimq=10):
        """Run trim_galore in single mode
        
        Returns
        -------
        tuple
            returns a tuple with two values. First is status of the command; True for success and False for failiure.
            second item in the tupple is the path to the qc corrected file
        """
        print ("Running bbduk single")
        
        #default output dir
        fnameSuffix="_bbduk"
        outDir=getFileDirectory(fastqFilePath)
        outFileName= getFileBaseName(fastqFilePath)+fnameSuffix+".fastq"
        outFilePath=os.path.join(outDir,outFileName)
        
        bbdukCmd=['bbduk.sh','-Xmx1g','in='+fastqFilePath,'out1='+outFilePath,'ref='+pathToAdapters,'ktrim='+ktrim,'k='+str(k),'mink='+str(mink),'hdist='+str(hdist),'qtrim='+qtrim,'trimq='+str(trimq),'threads='+str(proc)]
        print("Executing: "+" ".join(bbdukCmd))
        
        try:
            for output in executeCommand(bbdukCmd):
                print (output)
        except subprocess.CalledProcessError as e:
            print ("Error in command")
            return False,""
        
        
        #check if file exists
        if not checkFilesExists(outFilePath):
            print ("ERROR in running"+ self.programName)
            return False,""
        
        return True,outFilePath
        
        
            
    def runBBdukPaired(self,fastqFile1Path,fastqFile2Path,pathToAdapters="",proc="auto",ktrim='r',k=23,mink=11,hdist=1,qtrim='rl',trimq=10):
        """Function to run BBduk on paired data
        
        Returns
        -------
        tuple
            returns a tuple with three values. First is status of the command; True for success and False for failiure.
            second item in the tupple is the paths to the qc corrected file in order 1 and 2.
        """
        print ("Running bbduk paired")
        
        #default output dir
        fnameSuffix="_bbduk"
        outDir=getFileDirectory(fastqFile1Path)
        outFileName1=getFileBaseName(fastqFile1Path)+fnameSuffix+".fastq"
        outFileName2=getFileBaseName(fastqFile2Path)+fnameSuffix+".fastq"
        outFile1Path=os.path.join(outDir,outFileName1)
        outFile2Path=os.path.join(outDir,outFileName2)
        
        bbdukCmd=['bbduk.sh','-Xmx1g','in1='+fastqFile1Path,'in2='+fastqFile2Path,'out1='+outFile1Path,'out2='+outFile2Path,'ref='+pathToAdapters,'ktrim='+ktrim,'k='+str(k),'mink='+str(mink),'hdist='+str(hdist),'qtrim='+qtrim,'trimq='+str(trimq),'threads='+str(proc)]
        print("Executing: "+" ".join(bbdukCmd))
        
        
        try:
                for output in executeCommand(bbdukCmd):
                    print (output)
        except subprocess.CalledProcessError as e:
                print ("Error in command")
                return False,"",""
        print("Exiting...")
        
        #check if file exists
        if not checkFilesExists(outFile1Path,outFile2Path):
            print ("ERROR in running"+ self.programName)
            return False,"",""
        
        return True,outFile1Path,outFile2Path