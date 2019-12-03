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
        self.passedArgumentList=parseUnixStyleArgs(self.validArgsList,kwargs)
        print(self.passedArgumentList)
        
            
    def run(self,sraOb):
        """Run this class' program.
        The function run() is consistent for all QC classess.
        This function returns a tuple containing the status and the file paths to the qc-corrected fastq files (in order 1,2 for paired)
        """
        return self.runTrimGalore(sraOb)
            
            
    def runTrimGalore(self,sraOb):
        #get layout
        if sraOb.layout=='PAIRED':
            fq1=sraOb.localfastq1Path
            fq2=sraOb.localfastq2Path
            return self.runTrimGalorePaired(fq1,fq2)
        else:
            return self.runTrimGaloreSingle(sraOb.localfastqPath)
            
        
            
    def runTrimGaloreSingle(self,fastqFilePath):
        """Run trim_galore in single mode
        
        Returns
        -------
        tuple
            returns a tuple with two values. First is status of the command; True for success and False for failiure.
            second item in the tupple is the path to the qc corrected file
        """
        print("Running trim_galore unpaired")
        
        '''
        outfile will be named as <fastqFileName>_trimmed.fq
        this will be renamed to <fastqFileName><suffix>.fastq
        '''
        fnameSuffix="_trimGalore"
        outFileName=getFileBaseName(fastqFilePath)+"_trimmed.fq"
        newOutFileName=getFileBaseName(fastqFilePath)+fnameSuffix+".fastq"
        
        outDir=""
        trimGaloreCmd=['trim_galore']
        #check if out dir is specified
        if '-o' not in self.passedArgumentList:            
            #default output dir
            outDir=os.path.split(fastqFilePath)[0]
            trimGaloreCmd.extend(['-o',outDir])
        else:
            outDir=self.passedArgumentList[self.passedArgumentList.index('-o')+1]
        trimGaloreCmd.extend(self.passedArgumentList)
        trimGaloreCmd.extend([fastqFilePath])

        print("Executing: "+" ".join(trimGaloreCmd))
        try:
            for output in executeCommand(trimGaloreCmd):
                print (output)
        except subprocess.CalledProcessError as e:
            print ("Error in command")
            return False,""
        
        outFilePath=os.path.join(outDir,outFileName)
        newOutFilePath=os.path.join(outDir,newOutFileName)
        #rename the files
        rename_Cmd=['mv',outFilePath,newOutFilePath]
        print("Executing:"+ " ".join(rename_Cmd))
        if getCommandReturnValue(rename_Cmd)!=0:
            print("Error in moving files")
            return False,"",""
        #check if file exist
        if not checkFilesExists(newOutFilePath):
            print ("ERROR in running"+ self.programName)
            return False,""
        
        return True,newOutFilePath
        
            
    def runTrimGalorePaired(self,fastqFile1Path,fastqFile2Path):
        """Run trim_galore on paired data
        
        Returns
        -------
        tuple
            returns a tuple with three values. First is status of the command; True for success and False for failiure.
            second item in the tupple is the paths to the qc corrected file in order 1 and 2.
        """
        print ("Running trim_galore paired")
        
        '''
        out put files will be written as 
        <file name>_val_1.fq and <file name>_val_2.fq
        change this to <file name>_<suffix>.fastq and <file name>_<suffix>.fastq
        '''        
        fnameSuffix="_trimGalore"
        outFile1Name=getFileBaseName(fastqFile1Path)+"_val_1.fq"
        outFile2Name=getFileBaseName(fastqFile2Path)+"_val_2.fq"
        newOutFile1Name=getFileBaseName(fastqFile1Path)+fnameSuffix+".fastq"
        newOutFile2Name=getFileBaseName(fastqFile2Path)+fnameSuffix+".fastq"
        
        
        outDir="" #the out put directory
        trimGaloreCmd=['trim_galore']
        #check if out dir is specified
        if '-o' not in self.passedArgumentList:            
            #default output dir
            outDir=os.path.split(fastqFile1Path)[0]
            trimGaloreCmd.extend(['-o',outDir])
        else:
            outDir=self.passedArgumentList[self.passedArgumentList.index('-o')+1]
            
        trimGaloreCmd.extend(self.passedArgumentList)
        trimGaloreCmd.extend(['--paired',fastqFile1Path,fastqFile2Path])
        print("Executing: "+" ".join(trimGaloreCmd))
        try:
            for output in executeCommand(trimGaloreCmd):
                print (output)
        except subprocess.CalledProcessError as e:
            print ("Error in command")
            return False,"",""
        
        outFile1Path=os.path.join(outDir,outFile1Name)
        outFile2Path=os.path.join(outDir,outFile2Name)
        newOutFile1Path=os.path.join(outDir,newOutFile1Name)
        newOutFile2Path=os.path.join(outDir,newOutFile2Name)
        #rename the files
        rename_Cmd=['mv',outFile1Path,newOutFile1Path]
        print("Executing:"+ " ".join(rename_Cmd))
        if getCommandReturnValue(rename_Cmd)!=0:
            print("Error in moving files")
            return False,"",""
        rename_Cmd=['mv',outFile2Path,newOutFile2Path]
        if getCommandReturnValue(rename_Cmd)!=0:
            print("Error in moving files")
            return False,"",""
        
        #check if files exist
        if not checkFilesExists(newOutFile1Path,newOutFile2Path):
            print ("ERROR in running"+ self.programName)
            return False,"",""
            
        return True,newOutFile1Path,newOutFile2Path
            
            

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