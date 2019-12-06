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
    def __init__(self,**kwargs):
        """
        Parameters
        ----------
        kwargs:
            bbduk.sh arguments.
        """
        #run super to inherit parent class properties
        super().__init__() 
        self.programName="bbduk.sh"
        self.depList=[self.programName]
        #check if program exists
        if not checkDep(self.depList):
            raise Exception("ERROR: "+ self.programName+" not found.")
            
        
        self.validArgsList=['in','in2','ref','literal','touppercase','interleaved','qin','reads','copyundefined',
                            'samplerate','samref','out','out2','outm','outm2','outs','stats','refstats','rpkm',
                            'dump','duk','nzo','overwrite','showspeed','ziplevel','fastawrap','qout','statscolumns',
                            'rename','refnames','trd','ordered','maxbasesout','maxbasesoutm','','json','bhist','qhist',
                            'qchist','aqhist','bqhist','lhist','phist','gchist','ihist','gcbins','maxhistlen','histbefore',
                            'ehist','qahist','indelhist','mhist','idhist','idbins','varfile','vcf','ignorevcfindels',
                            'k','rcomp','maskmiddle','minkmerhits','minkmerfraction','mincovfraction','hammingdistance',
                            'qhdist','editdistance','hammingdistance2','qhdist2','editdistance2','forbidn','removeifeitherbad',
                            'trimfailures','findbestmatch','skipr1','skipr2','ecco','recalibrate','sam','le.','amino',
                            'threads','prealloc','monitor','minrskip','maxrskip','rskip','qskip','speed','ktrim','kmask',
                            'maskfullycovered','ksplit','mink','qtrim','trimq','trimclip','minlength','mlf','maxlength',
                            'minavgquality','maqb','minbasequality','maxns','mcb','ottm','tp','tbo','strictoverlap',
                            'minoverlap','mininsert','tpe','forcetrimleft','forcetrimright','forcetrimright2',
                            'forcetrimmod','restrictleft','restrictright','mingc','maxgc','gcpairs','tossjunk',
                            'swift','chastityfilter','barcodefilter','barcodes','xmin','ymin','xmax','ymax','trimpolya',
                            'trimpolygleft','trimpolygright','trimpolyg','filterpolyg','pratio','plen','entropy','entropywindow',
                            'entropyk','minbasefrequency','entropytrim','entropymask','entropymark','cardinality',
                            'cardinalityout','loglogk','loglogbuckets','-Xmx','-eoom','-da']
        
        self.passedArgumentDict=kwargs
            
            
            
    def performQC(self,sraOb,outFileSuffix="_bbduk",overwrite=True,**kwargs):
        """Run bbduk on fastq files specified by the sraOb
        
        Parameters
        ----------
        arg1: SRA
            an SRA object
        arg2: string
            Suffix for output file name
        arg3: bool
            overwrite existing files
        arg3: dict
            options passed to bbduk
            
        Returns
        tuple
            Returns the path of fastq files after QC. tuple has one item for single end files and 2 for paired.
        """
        if sraOb.layout=='PAIRED':
            fq1=sraOb.localfastq1Path
            fq2=sraOb.localfastq2Path
            #append input and output options
            outDir=sraOb.location
            outFileName1=getFileBaseName(fq1)+outFileSuffix+".fastq"
            outFileName2=getFileBaseName(fq2)+outFileSuffix+".fastq"
            outFile1Path=os.path.join(outDir,outFileName1)
            outFile2Path=os.path.join(outDir,outFileName2)
            
            newOpts={"in":fq1,"in2":fq2,"out":outFile1Path,"out2":outFile2Path}
            mergedOpts={**kwargs,**newOpts}
            
            #run bbduk
            if self.runBBduk(**mergedOpts):
                if checkFilesExists(outFile1Path,outFile2Path):
                    return(outFile1Path,outFile2Path)
            return("",)
            
            
        else:
            fq=sraOb.localfastqPath
            #append input and output options
            outDir=sraOb.location
            outFileName=getFileBaseName(fq)+outFileSuffix+".fastq"
            outFilePath=os.path.join(outDir,outFileName)
            newOpts={"in":fq,"out":outFilePath}
            mergedOpts={**kwargs,**newOpts}
            
            #run bbduk
            if self.runBBduk(**mergedOpts):
                if checkFilesExists(outFilePath):
                    return(outFilePath,)
            return("",)
                    
    
    
    
    
    def runBBduk(self,**kwargs):
        """Wrapper to run bbduk.sh
        """
        #override existing arguments
        mergedArgsDict={**self.passedArgumentDict,**kwargs}
        
        #create command to run
        bbduk_Cmd=["bbduk.sh"]
        
        #bbduk.sh follows java style arguments
        bbduk_Cmd.extend(parseJavaStyleArgs(self.validArgsList,mergedArgsDict))
        print("Executing:"+" ".join(bbduk_Cmd))
        
        #start ececution
        log=""
        try:
            for output in executeCommand(bbduk_Cmd):
                #print (output)    
                log=log+str(output)
            #save to a log file
        except subprocess.CalledProcessError as e:
            print ("Error in command...\n"+str(e))
            #save error to error.log file
            return False        
        #return status
        return True
    
    
    """
    index 
    bbsplit.sh ref_x=/home/usingh/work/urmi/hoap/test/bowtieIndex/euk_combined_rRNA.fa path=./xds
    
    clean
    bbsplit.sh in1=reads1.fq in2=reads2.fq ref=path_to_ref outu1=clean1.fq outu2=clean2.fq
    """
    
    def performCleaning(self,sraOb,bbsplitIndex,outFileSuffix="_bbsplit",overwrite=True,**kwargs):
        """
        Remove contaminated reads mapping to given reference using bbsplit
        
        Parameters
        ----------
        arg1: SRA
            an SRA object
        arg2: string
            Path to bbsplit index or fasta file which will generate index
        arg3: string
            Suffix for output file name
        arg4: bool
            overwrite existing files
        arg5: dict
            options passed to bbsplit
            
        Returns
        tuple
            Returns the path of fastq files after QC. tuple has one item for single end files and 2 for paired.
        """
        
        #check index
        indexPath=""
        if not checkPathsExists(bbsplitIndex):
            #index folder doesn't exist
            #check if input is path to fasta
            if not checkFilesExists(bbsplitIndex):
                print("Error: Please check bbspli index")
                return ("",)
            #check if index folder "ref" exists in this directory
            indexPath=os.path.join(getFileDirectory(bbsplitIndex),"ref")
            if checkPathsExists(indexPath):
                print("Using bbsplit index: "+indexPath)
            else:
                #create new index
                print("Creating new index"+indexPath)
                newOpts={"ref":bbsplitIndex,"path":getFileDirectory(bbsplitIndex)}
                mergedOpts={**kwargs,**newOpts}
                #run bbduk
                if not self.runBBsplit(**mergedOpts):
                    print("Error creating bbsplit index.")
                    return ("",)
                if not checkPathsExists(indexPath):
                    print("Error creating bbsplit index.")
                    return ("",)
        else:
            indexPath=bbsplitIndex
                
        
                
        
        
        if sraOb.layout=='PAIRED':
            fq1=sraOb.localfastq1Path
            fq2=sraOb.localfastq2Path
            #append input and output options
            outDir=sraOb.location
            outFileName1=getFileBaseName(fq1)+outFileSuffix+".fastq"
            outFileName2=getFileBaseName(fq2)+outFileSuffix+".fastq"
            outFile1Path=os.path.join(outDir,outFileName1)
            outFile2Path=os.path.join(outDir,outFileName2)
            
            newOpts={"in":fq1,"in2":fq2,"outu1":outFile1Path,"outu2":outFile2Path,"ref":indexPath}
            mergedOpts={**kwargs,**newOpts}
            
            #run bbduk
            if self.runBBsplit(**mergedOpts):
                if checkFilesExists(outFile1Path,outFile2Path):
                    return(outFile1Path,outFile2Path)
            return("",)
            
            
        else:
            fq=sraOb.localfastqPath
            #append input and output options
            outDir=sraOb.location
            outFileName=getFileBaseName(fq)+outFileSuffix+".fastq"
            outFilePath=os.path.join(outDir,outFileName)
            newOpts={"in":fq,"outu":outFilePath,"ref":indexPath}
            mergedOpts={**kwargs,**newOpts}
            
            #run bbduk
            if self.runBBsplit(**mergedOpts):
                if checkFilesExists(outFilePath):
                    return(outFilePath,)
            
            return("",)
    
    
    
    def runBBsplit(self,**kwargs):
        """wrapper to run bbsplit
        """
        #override existing arguments
        mergedArgsDict={**self.passedArgumentDict,**kwargs}
        
        #create command to run
        bbduk_Cmd=["bbsplit.sh"]
        
        #bbduk.sh follows java style arguments
        bbduk_Cmd.extend(parseJavaStyleArgs(self.validArgsList,mergedArgsDict))
        print("CCCCCDDDDDDDD Executing:"+" ".join(bbduk_Cmd))
        
        #start ececution
        log=""
        try:
            for output in executeCommand(bbduk_Cmd):
                #print (output)    
                log=log+str(output)
            #save to a log file
        except subprocess.CalledProcessError as e:
            print ("Error in command...\n"+str(e))
            #save error to error.log file
            return False        
        #return status
        return True
    
    
    
    
    
    
        