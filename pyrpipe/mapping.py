#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 19:53:42 2019

@author: usingh
contains classes of RNA-Seq mapping programs
"""

from pyrpipe.pyrpipe_utils import *
from pyrpipe.pyrpipe_engine import *

class Aligner:
    def __init__(self):
        self.category="Alignement"
        self.passedArgumentDict={}
    
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
        #check if hisat2 exists
        if not checkDep([self.programName]):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        self.validArgsList=['-x','-1','-2','-U','--sra-acc','-S','-q','--qseq','-f','-r','-c','-s',
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
        
        
        #initialize the passed arguments
        self.passedArgumentDict=kwargs
        
        #if index is passed, update the passed arguments
        if len(hisat2Index)>0 and checkHisatIndex(hisat2Index):
            print("HISAT2 index is: "+hisat2Index)
            self.hisat2Index=hisat2Index
            self.passedArgumentDict['-x']=self.hisat2Index
        else:
            print("No Hisat2 index provided. Please build index now to generate an index using buildHisat2Index()....")
            
        
        
            
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
        status=executeCommand(hisat2Build_Cmd)
        if not status:
            printBoldRed("hisatBuild failed")
            return False
        
        #check if sam file is present in the location directory of sraOb
        if not checkHisatIndex(os.path.join(indexPath,indexName)):
            printBoldRed("hisatBuild failed")
            return False
        
        #set the index path
        self.hisat2Index=os.path.join(indexPath,indexName)
        self.passedArgumentDict['-x']=self.hisat2Index
        
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
        
        
        #create path to output sam file
        outSamFile=os.path.join(sraOb.location,sraOb.srrAccession+outSamSuffix+".sam")
        
        """
        Handle overwrite
        """
        overwrite=True
        if not overwrite:
            #check if file exists. return if yes
            if os.path.isfile(outSamFile):
                print("The file "+outSamFile+" already exists. Exiting..")
                return outSamFile
        
        #find layout and fq file paths
        if sraOb.layout == 'PAIRED':
            newOpts={"-1":sraOb.localfastq1Path,"-2":sraOb.localfastq2Path,"-S":outSamFile}
        else:
            newOpts={"-U":sraOb.localfastqPath,"-S":outSamFile}
        
        #add input files to kwargs, overwrite kwargs with newOpts
        mergedOpts={**kwargs,**newOpts}
        
        #call runHisat2
        status=self.runHisat2(**mergedOpts)
        
        if status:
            #check if sam file is present in the location directory of sraOb
            if checkFilesExists(outSamFile):
                return outSamFile
        else:
            return ""
            
        
    def runHisat2(self,**kwargs):
        """Wrapper for running hisat2.
        Run HISAT2 using and SRA object and produce .bam file as result. The HISAT2 index used will be self.hisat2Index.
        All output will be written to SRA.location by default.
        
        Parameters
        ----------
        arg1: dict
            arguments to pass to hisat2. This will override parametrs already existing in the self.passedArgumentList list but NOT replace them.
            
        Returns
        -------
        bool:
                Returns the status of hisat2. True is passed, False if failed.
        """
        
        #check for a valid index
        if not self.checkHisat2Index():
            raise Exception("ERROR: Invalid HISAT2 index. Please run build index to generate an index.")
            
        #override existing arguments
        mergedArgsDict={**self.passedArgumentDict,**kwargs}
       
        hisat2_Cmd=['hisat2']
        #add options
        hisat2_Cmd.extend(parseUnixStyleArgs(self.validArgsList,mergedArgsDict))        
        
        #execute command
        cmdStatus=executeCommand(hisat2_Cmd)
        if not cmdStatus:
            print("hisat2 failed:"+" ".join(hisat2_Cmd))
     
        #return status
        return cmdStatus
        
        
    
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
    def __init__(self,bowtie2Index,**kwargs):
        """Bowtie2 constructor. Initialize star's index and other parameters.
        """       
        
        super().__init__() 
        self.programName="bowtie2"
        self.depList=[self.programName]        
        if not checkDep(self.depList):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        self.validArgsList=['-x','-1','-2','-U','--interleaved','-S','-b','-q','--tab5','--tab6','--qseq','-f','-r','-F','-c','-s','-u','-5','-3',
                            '--trim-to','--phred33','--phred64','--int-quals','--very-fast','--fast',
                            '--sensitive','--very-sensitive','--very-fast-local','--fast-local',
                            '--sensitive-local','--very-sensitive-local','-N','-L','-i','--n-ceil',
                            '--dpad','--gbar','--ignore-quals','--nofw','--norc','--no-1mm-upfront',
                            '--end-to-end','--local','--ma','--mp','--np','--rdg','--rfg','--score-min',
                            '-k','-a','-D','-R','-I','-X','--fr','--rf','--ff','--no-mixed','--no-discordant',
                            '--dovetail','--no-contain','--no-overlap','--align-paired-reads','--preserve-tags',
                            '-t','--un','--al','--un-conc','--al-conc','--un-gz','--quiet','--met-file',
                            '--met-stderr','--met','--no-unal','--no-head','--no-sq','--rg-id','--rg',
                            '--omit-sec-seq','--sam-no-qname-trunc','--xeq','--soft-clipped-unmapped-tlen',
                            '-p','--threads','--reorder','--mm','--qc-filter','--seed','--non-deterministic',
                            '--version','-h','--help']
        
        #initialize the passed arguments
        self.passedArgumentDict=kwargs
        
        #if index is passed, update the passed arguments
        if len(bowtie2Index)>0 and checkBowtie2Index(bowtie2Index):
            print("Bowtie2 index is: "+bowtie2Index)
            self.bowtie2Index=bowtie2Index
            self.passedArgumentDict['-x']=self.bowtie2Index
        else:
            print("No Bowtie2 index provided. Please build index now to generate an index...")
        
        
    
    
    def performAlignment(self,sraOb,outSamSuffix="_bt2",overwrite=True,**kwargs):
        """Function to perform alignment using self object and the provided sraOb.
        
        Parameters
        ----------
        arg1: SRA object
            An object of type SRA. The path to fastq files will be obtained from this object.
        arg2: string
            Suffix for the output sam file
        arg3: dict
            Options to pass to bowtie.
        """
        
        #create path to output sam file
        outFile=os.path.join(sraOb.location,sraOb.srrAccession+outSamSuffix+".sam")
                    
        """
        Handle overwrite
        """
        overwrite=True
        if not overwrite:
            #check if file exists. return if yes
            if os.path.isfile(outFile):
                print("The file "+outFile+" already exists. Exiting..")
                return outFile
        
        #find layout and fq file paths
        if sraOb.layout == 'PAIRED':
            newOpts={"-1":sraOb.localfastq1Path,"-2":sraOb.localfastq2Path,"-S":outFile}
        else:
            newOpts={"-U":sraOb.localfastqPath,"-S":outFile}
        
        #add input files to kwargs, overwrite kwargs with newOpts
        mergedOpts={**kwargs,**newOpts}
        
        status=self.runBowTie2(**mergedOpts)
        
        if status:
            #check if sam file is present in the location directory of sraOb
            if checkFilesExists(outFile):
                return outFile
        else:
            return ""
        
        
        
    
    def runBowTie2(self,**kwargs):
        """Wrapper for running bowtie2.
        
        ----------
        arg1: dict
            arguments to pass to bowtie2. This will override parametrs already existing in the self.passedArgumentList list but NOT replace them.
            
        Returns
        -------
        bool:
                Returns the status of bowtie2. True is passed, False if failed.
        """
        
        #check for a valid index
        if not self.checkIndex():
            raise Exception("ERROR: Invalid Bowtie2 index. Please run build index to generate an index.")
        
        #override existing arguments
        mergedArgsDict={**self.passedArgumentDict,**kwargs}
            
        bowtie2_Cmd=['bowtie2']
        bowtie2_Cmd.extend(parseUnixStyleArgs(self.validArgsList,mergedArgsDict))
        
        print("Executing:"+" ".join(bowtie2_Cmd))
        
        #start ececution
        status=executeCommand(bowtie2_Cmd)
        if not status:
            printBoldRed("bowtie2 failed")
        return status
    
    
    def checkIndex(self):
        if hasattr(self,'bowtie2Index'):
            return(checkBowtie2Index(self.bowtie2Index))
        else:
            return False
        
    
    
    
