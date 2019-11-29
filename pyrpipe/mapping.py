#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 19:53:42 2019

@author: usingh
contains classes of RNA-Seq mapping programs
"""

from pyrpipe.myutils import *

class Aligner:
    def __init__(self):
        self.category="Alignement"

class Hisat2:
    def __init__(self,hisat2Index="",*args):
        """HISAT2 constructor. Initialize hisat2's index and other parameters.
        Parameters
        ----------
        
        """ 
        super().__init__() 
        self.programName="hisat2"
        #check if hisat2 exists
        if not checkDep([self.programName]):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        
        #check index exists
        if len(hisat2Index)>0 and checkHisatIndex(hisat2Index):
            print("Found HISAT2 index files.")
            self.hisat2Index=hisat2Index
        elif args :
            print("Building HISAT2 index with "+ str(args))
            indexPath,indexName=os.path.split(hisat2Index)
            self.buildHisat2Index(indexPath,indexName,*args)
        else:
            raise Exception("ERROR: Invalid HISAT2 index. Please run build index to generate an index.")
            
    def buildHisat2Index(self,indexPath,indexName,*args):
        print("Building hisat index...")
        
        
    
    def runHisat2(self,sraOb,outSamSuffix="_hisat2",**kwargs):
        """Run HISAT2 using and SRA object and produce .bam file as result. The HISAT2 index used will be self.hisat2Index.
        All output will be written to SRA.location by default.
        
        Parameters
        ----------
        arg1: SRA object
            An object of type SRA which contains valid fastq files
        
        arg2: dict
            arguments to be passed to hisat2
        """
        
        #check for a valid index
        if not self.checkHisat2Index():
            raise Exception("ERROR: Invalid HISAT2 index. Please run build index to generate an index.")
            
         #save information about the SRAobject
         #self.SRAob=SRAob
         
        #scan for prefetch arguments
        hisatArgsList=['-p','--dta-cufflinks']
        
        outSamFile=os.path.join(sraOb.location,sraOb.srrAccession+outSamSuffix+".sam")
        
        
        #check if file exists. return if yes
        if os.path.isfile(outSamFile):
            print("The file "+outSamFile+" already exists. Exiting..")
            return True,outSamFile
            
        hisat2_Cmd=['hisat2']
        hisat2_Cmd.extend(parseUnixStyleArgs(hisatArgsList,kwargs))
        hisat2_Cmd.extend(['-x',self.hisat2Index])
        if sraOb.layout == 'PAIRED':
            hisat2_Cmd.extend(['-1',sraOb.localfastq1Path])
            hisat2_Cmd.extend(['-2',sraOb.localfastq2Path])
        else:
            hisat2_Cmd.extend(['-U',sraOb.localfastqPath])
        #save output to the sraob location folder 
        hisat2_Cmd.extend(['-S',outSamFile])
        print("Executing:"+" ".join(hisat2_Cmd))
        
        #start ececution
        log=""
        try:
            for output in executeCommand(hisat2_Cmd):
                print (output)    
                log=log+str(output)
            #save to a log file

        except subprocess.CalledProcessError as e:
            print ("Error in command...\n"+str(e))
            #save error to error.log file
            return False,"NA"
        
        #check if sam file is present in the location directory of sraOb
        if not checkFilesExists(outSamFile):
            return False,"NA"
        
        #return status and path to output sam
        return True,outSamFile
        
        
    
    def checkHisat2Index(self):
        if hasattr(self,'hisat2Index'):
            return(checkHisatIndex(self.hisat2Index))
        else:
            return False





class Star:
    def __init__(self,starIndex):
        """STAR constructor. Initialize star's index and other parameters.
        """
        
        self.programName="star"
        #check if hisat2 exists
        if not checkDep([self.programName]):
            raise Exception("ERROR: "+ self.programName+" not found.")
            

    
class Samtools:
    def __init__(self):
        self.programName="samtools"
        #check if hisat2 exists
        if not checkDep([self.programName]):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
    def samToBam(self,samFile,proc,deleteSam=False):
        """Convert sam file to a bam file. Output bam file will have same name as input sam.
        """
        fname=samFile.split('.sam')[0]
        outBamFile=fname+'.bam'
        samToBamCmd=['samtools','view','-@',str(proc),'-b','-o',outBamFile,samFile]
        print("Executing: "+" ".join(samToBamCmd))
        try:
            for output in executeCommand(samToBamCmd):
                print(output)
        except subprocess.CalledProcessError as e:
            print ("Error in command")
            return -1
        print("sam2bam finished")
    
        #deletesamfile
        if deleteSam:
            delSamCmd=['rm',samFile]
            print("Deleting sam file...")
            try:
                    for output in executeCommand(delSamCmd):
                            print (output)
            except subprocess.CalledProcessError as e:
                print ("Error in command")
                return -1
        
        #check if bam file exists
        if not checkFilesExists(outBamFile):
            return False,"NA"
        #return status and file
        return True,outBamFile
        
        
        
        
    #sort bam file.output will be bamFile_sorted.bam
    def sortBam(self,bamFile,proc,deleteOriginalBam=False):
        fname=bamFile.split('.bam')[0]
        outSortedBamFile=fname+"_sorted.bam"
        bamSortCmd=['samtools','sort','-o',outSortedBamFile,'-@',str(proc),bamFile]
        print("Executing: "+" ".join(bamSortCmd))
        try:
            for output in executeCommand(bamSortCmd):
                print (output)
        except subprocess.CalledProcessError as e:
            print ("Error in command")
            return False
        print("SamSort finished")
    
        if deleteOriginalBam:
            delBamCmd=['rm',bamFile]
            print("Deleting unsorted bam file...")
            try:
                for output in executeCommand(delBamCmd):
                         print (output)
            except subprocess.CalledProcessError as e:
                 print ("Error in command")
                 return False
        
        #check if bam file exists
        if not checkFilesExists(outSortedBamFile):
            return False,"NA"
        #return status and file
        return True,outSortedBamFile
    
    def samToSortedBam(self,samFile,proc,deleteSam=False,deleteOriginalBam=False):
        """Convert sam file to bam and sort the bam file.
        """
        
        sam2BamStatus=self.samToBam(samFile,proc,deleteSam)
        
        if not sam2BamStatus[0]:
            raise Exception("ERROR in converting sam to bam: "+samFile)
        
        bamSortStatus=self.sortBam(sam2BamStatus[1],proc,deleteOriginalBam)
        
        if not bamSortStatus[0]:
            raise Exception("ERROR sorting bam file: "+sam2BamStatus[1])
        
        return bamSortStatus[1]

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        