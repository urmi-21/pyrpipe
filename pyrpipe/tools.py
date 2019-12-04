#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 14:54:22 2019

@author: usingh
"""
from pyrpipe.myutils import *

class RNASeqTools:
    def __init__(self):
        self.category="RNASeqTools"
    
    
class Samtools(RNASeqTools):
    def __init__(self,**kwargs):
        self.programName="samtools"
        #check if hisat2 exists
        if not checkDep([self.programName]):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        self.validArgsList=['-b','-C','-1','-u','-h','-H','-c','-o','-U','-t','-L','-r',
                            '-R','-q','-l','-m','-f','-F','-G','-s','-M','-x','-B','-?','-S','-O','-T','-@']
        
        self.passedArgumentDict=kwargs
        
        
        
    def samToBam(self,samFile,outFileSuffix="",deleteSam=False,**kwargs):
        """Convert sam file to a bam file. 
        Output bam file will have same name as input sam.
        
        Returns
        -------
        string
                Returns the path to the bam file. Returns empty string if operation failed.
        """        
        
        outDir=getFileDirectory(samFile)
        fname=getFileBaseName(samFile)
        #output will be outBamFile
        outBamFile=os.path.join(outDir,fname+outFileSuffix+'.bam')
        
        newOpts={"--":(samFile,),"-o":outBamFile,"-b":""}
        mergedOpts={**kwargs,**newOpts}
        
        status=self.runSamtools("view",**mergedOpts)
                
        if not status:
            print("Sam to bam failed for:"+samFile)
            return ""
        
        #check if bam file exists
        if not checkFilesExists(outBamFile):
            return ""
        
        #deletesamfile
        if deleteSam:
            if not deleteFileFromDisk(samFile):
                print("Error deleting sam file:"+samFile)
                
        #return path to file
        return outBamFile
        
        
        
        
    #sort bam file.output will be bamFile_sorted.bam
    def sortBam(self,bamFile,outFileSuffix="",deleteOriginalBam=False,**kwargs):
        """Sorts an input bam file. Outpufile will end in _sorted.bam
        
        Returns
        -------
        string
                Returns path to the sorted bam file. Returns empty string if operation failed.
        
        """
        
        outDir=getFileDirectory(bamFile)
        fname=getFileBaseName(bamFile)
        #output will be outBamFile
        outSortedBamFile=os.path.join(outDir,fname+outFileSuffix+'_sorted.bam')
        
        newOpts={"--":(bamFile,),"-o":outSortedBamFile}
        mergedOpts={**kwargs,**newOpts}
        
        status=self.runSamtools("sort",**mergedOpts)
        
        if not status:
            print("Bam sort failed for:"+bamFile)
            return ""
        
        #check if bam file exists
        if not checkFilesExists(outSortedBamFile):
            return ""

        if deleteOriginalBam:
            if not deleteFileFromDisk(bamFile):
                print("Error deleting sam file:"+bamFile)
                
        #return path to file
        return outSortedBamFile
    
    def samToSortedBam(self,samFile,outFileSuffix="",deleteSam=False,deleteOriginalBam=False,**kwargs):
        """Convert sam file to bam and sort the bam file.
        
        Returns
        -------
        string
                Returns path to the sorted bam file. Returns empty string if operation failed.
        """
        
        sam2BamFile=self.samToBam(samFile,deleteSam=deleteSam,**kwargs)
        
        if not sam2BamFile:
            return ""
            

        bamSorted=self.sortBam(sam2BamFile,outFileSuffix,deleteOriginalBam,**kwargs)
        
        if not bamSorted:
            return ""
        
        return bamSorted
    
    
    def runSamtools(self,subCommand,**kwargs):
        """A wrapper to run samtools.
        
        Parameters
        ----------
        arg1: dict
            arguments to pass to samtools. This will override parametrs already existing in the self.passedArgumentDict list but NOT replace them.
            
        Returns
        -------
        bool:
                Returns the status of samtools. True is passed, False if failed.
        """
            
        #override existing arguments
        mergedArgsDict={**self.passedArgumentDict,**kwargs}
       
        samtools_Cmd=['samtools',subCommand]
        #add options
        samtools_Cmd.extend(parseUnixStyleArgs(self.validArgsList,mergedArgsDict))
                
        print("Executing:"+" ".join(samtools_Cmd))
        
        
        #start ececution
        log=""
        try:
            for output in executeCommand(samtools_Cmd):
                #print (output)    
                log=log+str(output)
            #save to a log file
            
        except subprocess.CalledProcessError as e:
            print ("Error in command...\n"+str(e))
            #save error to error.log file
            return False        
        #return status
        return True
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        