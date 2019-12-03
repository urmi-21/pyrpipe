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
    
    def performAlignment(self):
        pass

class Hisat2:
    def __init__(self,hisat2Index="",**kwargs):
        """HISAT2 constructor. Initialize hisat2's index and other parameters.
        Parameters
        ----------
        hisat2Index string
            path to q histat2 index (note -x is ommited from validArgsList). This index will be used when hisat is invoked.
        dict
            parameters passed to the hisat2 program. These parameters could be overridden later when running hisat.
        ----------
        
        """ 
        super().__init__() 
        self.programName="hisat2"
        
        self.validArgsList=['-1','-2','-U','--sra-acc','-S','-q','--qseq','-f','-r','-c','-s',
                            '-u','-5','-3','--phred33','--phred64','--int-quals',
                            '--sra-acc','--n-ceil','--ignore-quals','--nofw','--norc','--pen-cansplice',
                            '--pen-noncansplice','--pen-canintronlen','--pen-noncanintronlen','--min-intronlen'
                            ,'--max-intronlen','--known-splicesite-infile','--novel-splicesite-outfile',
                            '--novel-splicesite-infile','--no-temp-splicesite','--no-spliced-alignment',
                            '--rna-strandness','--tmo','--dta','--dta-cufflinks','--avoid-pseudogene',
                            '--no-templatelen-adjustment','--mp','--sp','--no-softclip','--np','--rdg',
                            '--rfg','--score-min','-k','-I','-X','--fr','--rf','--ff','--no-mixed',
                            '--no-discordant','-t','--un','--al','--un-conc','--al-conc','--un-gz',
                            '--summary-file','--new-summary','--quiet','--met-file','--met-stderr',
                            '--met','--no-head','--no-sq','--rg-id','--rgit-sec-seq','-o','-p',
                            '--reorder','--mm','--qc-filter','--seed','--non-deterministic',
                            '--remove-chrname','--add-chrname','--version']
        
        
        
        #check if hisat2 exists
        if not checkDep([self.programName]):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        
        #check index exists
        if len(hisat2Index)>0 and checkHisatIndex(hisat2Index):
            print("HISAT2 index is: "+hisat2Index)
            self.hisat2Index=hisat2Index
        else:
            print("No Hisat2 index provided. Please run build index now to generate an index....")
            
        #initialize the passed arguments
        self.passedArgumentDict=kwargs
        self.passedArgumentList=parseUnixStyleArgs(self.validArgsList,kwargs)
            
    def buildHisat2Index(self,indexPath,indexName,*args,**kwargs):
        """Build a hisat index with given parameters and saves the new index to self.hisat2Index.
        Parameters
        ----------
        arg1: string
            Path where the index will be created
        arg2: string
            A name for the index
        arg3: tuple
            Path to reference input files
        arg4: dict
            Parameters for the hisat2-build command
        
        Returns
        -------
        bool:
            Returns the status of hisat2-build
        """
        overwrite=True
        print("Building hisat index...")
        
        hisat2BuildValidArgsList=['-c','--large-index','-a','-p','--bmax','--bmaxdivn','--dcv','--nodc','-r','-3','-o',
                                  '-t','--localoffrate','--localftabchars','--snp','--haplotype','--ss','--exon',
                                  '--seed','-q','-h','--usage','--version']
        #create the out dir
        if not checkPathsExists(indexPath):
            if not mkdir(indexPath):
                print("ERROR in building hisat2 index. Failed to create index directory.")
                return False
        
        if not overwrite:
            #check if files exists
            if checkHisatIndex(os.path.join(indexPath,indexName)):
                print("Hisat2 index with same name already exists. Exiting...")
                return False
        
        hisat2Build_Cmd=['hisat2-build']
        #add options
        hisat2Build_Cmd.extend(parseUnixStyleArgs(hisat2BuildValidArgsList,kwargs))
        #add input files
        hisat2Build_Cmd.append(str(",".join(args)))
        #add dir/basenae
        hisat2Build_Cmd.append(os.path.join(indexPath,indexName))
        print("Executing:"+str(" ".join(hisat2Build_Cmd)))
        
        #start ececution
        log=""
        try:
            for output in executeCommand(hisat2Build_Cmd):
                print (output)    
                log=log+str(output)
            #save to a log file

        except subprocess.CalledProcessError as e:
            print ("Error in command...\n"+str(e))
            #save error to error.log file
            return ""
        
        #check if sam file is present in the location directory of sraOb
        if not checkHisatIndex(os.path.join(indexPath,indexName)):
            print("ERROR in building hisat2 index.")
            return False
        
        #set the index path
        self.hisat2Index=os.path.join(indexPath,indexName)
        
        #return the path to output sam
        return True
        
        
    def performAlignment(self,sraOb,outSamSuffix="_hisat2",**kwargs):
        """Function to perform alignment using self object and the provided sraOb.
        
        Parameters
        ----------
        arg1: SRA object
            An object of type SRA. The path to fastq files will be obtained from this object.
        arg2: string
            Suffix for the output sam file
        arg3: dict
            Options to pass to hisat2.
        """
        
        #find layout and fq file paths
        if sraOb.layout == 'PAIRED':
            pairedFlag=True
            fastqFileList=[sraOb.localfastq1Path,sraOb.localfastq2Path]
        else:
            pairedFlag=False
            fastqFileList=[sraOb.localfastqPath]
        
        #create path to output sam file
        outSamFile=os.path.join(sraOb.location,sraOb.srrAccession+outSamSuffix+".sam")
        
        
        #call runHisat2
        return self.runHisat2(fastqFileList,pairedFlag,outSamFile,**kwargs)
            
        
        
        
        
    
    #def runHisat2(self,sraOb,outSamSuffix="_hisat2",**kwargs):
    def runHisat2(self,fastqFileList,pairedFlag,outSamFile,**kwargs):
        """Run HISAT2 using and SRA object and produce .bam file as result. The HISAT2 index used will be self.hisat2Index.
        All output will be written to SRA.location by default.
        
        Parameters
        ----------
        arg1: list
            A list containing path to fastq files
        
        arg2: bool
            bool indicating whether data is paired. True --> Paired; False --> single
        
        arg3: string
            path to output sam file
        
        arg4: dict
            arguments to pass to hisat2. This will override parametrs already existing in the self.passedArgumentList list but NOT replace them.
            
        Returns
        -------
        string:
                Returns the path to the sam file created. If hisat2 fails, then returns the empty string.
        """
        
        #check for a valid index
        if not self.checkHisat2Index():
            raise Exception("ERROR: Invalid HISAT2 index. Please run build index to generate an index.")
        
        """
        Handle overwrite
        """
        overwrite=True
        if not overwrite:
            #check if file exists. return if yes
            if os.path.isfile(outSamFile):
                print("The file "+outSamFile+" already exists. Exiting..")
                return outSamFile
            
        #override existing arguments
        mergedArgsDict={**self.passedArgumentDict,**kwargs}
       
        
        hisat2_Cmd=['hisat2']
        #add options
        hisat2_Cmd.extend(parseUnixStyleArgs(self.validArgsList,mergedArgsDict))
        hisat2_Cmd.extend(['-x',self.hisat2Index])
        if pairedFlag:
            hisat2_Cmd.extend(['-1',fastqFileList[0]])
            hisat2_Cmd.extend(['-2',fastqFileList[1]])
        else:
            hisat2_Cmd.extend(['-U',fastqFileList[0]])
        #save output to the sraob location folder 
        hisat2_Cmd.extend(['-S',outSamFile])
        print("Executing:"+" ".join(hisat2_Cmd))
        
        
        
        #start ececution
        log=""
        try:
            for output in executeCommand(hisat2_Cmd):
                #print (output)    
                log=log+str(output)
            #save to a log file
            
        except subprocess.CalledProcessError as e:
            print ("Error in command...\n"+str(e))
            #save error to error.log file
            return ""
        
        #check if sam file is present in the location directory of sraOb
        if not checkFilesExists(outSamFile):
            return ""
        
        #return the path to output sam
        return outSamFile
        
        
    
    def checkHisat2Index(self):
        if hasattr(self,'hisat2Index'):
            return(checkHisatIndex(self.hisat2Index))
        else:
            return False

    



class Star:
    def __init__(self,starIndex):
        """STAR constructor. Initialize star's index and other parameters.
        """
        
        
        super().__init__() 
        self.programName="star"
        
        self.validArgsList=[]
        self.depList=[self.programName]        
        #check if hisat2 exists
        if not checkDep(self.depList):
            raise Exception("ERROR: "+ self.programName+" not found.")
            
            
            
    def performAlignment(self,sraOb,outSamSuffix="_star",**kwargs):
        """Function to perform alignment using self object and the provided sraOb.
        
        Parameters
        ----------
        arg1: SRA object
            An object of type SRA. The path to fastq files will be obtained from this object.
        arg2: string
            Suffix for the output sam file
        arg3: dict
            Options to pass to hisat2.
        """
        
        pass
            

class Bowtie2(Aligner):
    def __init__(self,bowtie2Index):
        """Bowtie2 constructor. Initialize star's index and other parameters.
        """       
        
        super().__init__() 
        self.programName="bowtie2"
        
        self.validArgsList=[]
        self.depList=[self.programName]        
        #check if hisat2 exists
        if not checkDep(self.depList):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        
        self.bowtie2Index=bowtie2Index
    
    
    def performAlignment(self,sraOb,outSamSuffix="_star",**kwargs):
        """Function to perform alignment using self object and the provided sraOb.
        
        Parameters
        ----------
        arg1: SRA object
            An object of type SRA. The path to fastq files will be obtained from this object.
        arg2: string
            Suffix for the output sam file
        arg3: dict
            Options to pass to hisat2.
        """
        
        pass
    
    def runBowTie2(self,sraOb,**kwargs):
        """Function to run bowtie2
        """
        outDir=sraOb.location
        outSamFile=os.path.join(outDir,"bt2.sam")
        unmapFname=os.path.join(outDir,sraOb.srrAccession+"norRNA.fastq")
        
        bowtie2_Cmd=[self.programName]
        bowtie2_Cmd.extend(["-p","10","--norc","-x",self.bowtie2Index,"-S",outSamFile,"--un",unmapFname])
        
        print("Executing:"+" ".join(bowtie2_Cmd))
        
    
    
    
class Samtools:
    def __init__(self):
        self.programName="samtools"
        #check if hisat2 exists
        if not checkDep([self.programName]):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
    def samToBam(self,samFile,proc,deleteSam=False):
        """Convert sam file to a bam file. Output bam file will have same name as input sam.
        
        Returns
        -------
        string
                Returns the path to the bam file. Returns empty string if operation failed.
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
            return ""
        
    
        #deletesamfile
        if deleteSam:
            delSamCmd=['rm',samFile]
            print("Deleting sam file...")
            try:
                    for output in executeCommand(delSamCmd):
                            print (output)
            except subprocess.CalledProcessError as e:
                print ("Error deleting sam...")
        
        #check if bam file exists
        if not checkFilesExists(outBamFile):
            return ""
        #return path to file
        return outBamFile
        
        
        
        
    #sort bam file.output will be bamFile_sorted.bam
    def sortBam(self,bamFile,proc,deleteOriginalBam=False):
        """Sorts an input bam file.
        
        Returns
        -------
        string
                Returns path to the sorted bam file. Returns empty string if operation failed.
        
        """
        fname=bamFile.split('.bam')[0]
        outSortedBamFile=fname+"_sorted.bam"
        bamSortCmd=['samtools','sort','-o',outSortedBamFile,'-@',str(proc),bamFile]
        print("Executing: "+" ".join(bamSortCmd))
        try:
            for output in executeCommand(bamSortCmd):
                print (output)
        except subprocess.CalledProcessError as e:
            print ("Error in command")
            return ""
        
    
        if deleteOriginalBam:
            delBamCmd=['rm',bamFile]
            print("Deleting unsorted bam file...")
            try:
                for output in executeCommand(delBamCmd):
                         print (output)
            except subprocess.CalledProcessError as e:
                 print ("Error deleting unsorted bam file...")
                 
        
        #check if bam file exists
        if not checkFilesExists(outSortedBamFile):
            return ""
        #return path to file
        return outSortedBamFile
    
    def samToSortedBam(self,samFile,proc,deleteSam=False,deleteOriginalBam=False):
        """Convert sam file to bam and sort the bam file.
        
        Returns
        -------
        string
                Returns path to the sorted bam file. Returns empty string if operation failed.
        """
        
        sam2BamFile=self.samToBam(samFile,proc,deleteSam)
        
        if not sam2BamFile:
            return ""
            
        bamSorted=self.sortBam(sam2BamFile,proc,deleteOriginalBam)
        
        if not bamSorted:
            return ""
        
        return bamSorted

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        