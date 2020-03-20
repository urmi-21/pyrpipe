#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 17:48:00 2019

@author: usingh
"""

from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe
import os
import math

class RNASeqQC:
    """This is an abstract parent class for fastq quality control programs.
    """
    def __init__(self):
        self.category="RNASeqQC"
        
    def perform_qc(self):
        pass

class Trimgalore(RNASeqQC):
    """This class represents trimgalore
    
    Parameters
    ----------
        
        kwargs:
            trim_galore arguments.
    """
    def __init__(self,threads=None):
        """
        threads: int
            Num threads to use
        """
        
        #run super to inherit parent class properties
        super().__init__() 
        self.programName="trim_galore"
        self.dep_list=[self.programName,'cutadapt']
        
        """
        self.valid_args=['--cores','-v','-q','--phred33','--phred64','--fastqc','--fastqc_args','-a','-a2',
                            '--illumina','--nextera','--small_rna','--consider_already_trimmed',
                            '--max_length','--stringency','-e','--gzip','--dont_gzip','--length',
                            '--max_n','--trim-n','-o','--no_report_file','--suppress_warn',
                            '--clip_R1','--clip_R2','--three_prime_clip_R1','--three_prime_clip_R2',
                            '--2colour','--path_to_cutadapt','--basename','-j','--hardtrim5','--hardtrim3',
                            '--clock','--polyA','--rrbs','--non_directional','--keep','--paired','-t',
                            '--retain_unpaired','-r1','-r2']
        """
        #check if deps exists
        if not pe.check_dependencies(self.dep_list):
            raise Exception("ERROR: "+ self.programName+" not found.")
            
        #initialize the passed arguments
        if not threads:
            #trimgalore recommends max 8 threads
            threads=8
        self.threads=threads
        
        
            
    def perform_qc(self,sra_object,out_dir="",out_suffix="_trimgalore",threads=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Function to perform qc using trimgalore.
        The function perform_qc() is consistent for all QC classess.
        
        Parameters
        ----------
        
        sra_object: SRA
            An SRA object whose fastq files will be used
        out_dir: str
            Path to output directory
        out_suffix: string
            Suffix for the output sam file
        threads: int
            Num threads to use
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to trimgalore. This will override the existing options 

            :return: Returns the path of fastq files after QC. tuple has one item for single end files and two for paired.
            :rtype: tuple
        """
        
        if not out_dir:
            out_dir=sra_object.location
        else:
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
        
        if not threads:
            threads=self.threads
        
        #create new options based on parametrs
        newOpts={}
        #get layout
        if sra_object.layout=='PAIRED':
            fq1=sra_object.localfastq1Path
            fq2=sra_object.localfastq2Path
            out_file1=os.path.join(out_dir,pu.get_file_basename(fq1)+out_suffix+".fastq")
            out_file2=os.path.join(out_dir,pu.get_file_basename(fq2)+out_suffix+".fastq")
            newOpts={"--paired":"","--":(fq1,fq2),"-o":out_dir,"--cores":str(threads)}
            mergedOpts={**newOpts,**kwargs}
            #run trimgalore
            self.run_trimgalore(verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
            """
            running trim galore will create two files named <input>_val_1.fq and <input>_val_2.fq
            move these files to the specified out files
            """
            oldFile1=os.path.join(out_dir,pu.get_file_basename(fq1)+"_val_1.fq")
            oldFile2=os.path.join(out_dir,pu.get_file_basename(fq2)+"_val_2.fq")
            
            pe.move_file(oldFile1,out_file1,verbose=False)
            pe.move_file(oldFile2,out_file2,verbose=False)
            
            if not pu.check_files_exist(out_file1,out_file2):
                print("Trimgalore failed")
                return ("",)
            return out_file1,out_file2
            
        else:
            fq=sra_object.localfastqPath
            out_file=os.path.join(out_dir, pu.get_file_basename(fq)+out_suffix+".fastq")
            #giving input arguments as a tuple "--":(fq,)
            newOpts={"--":(fq,),"-o":out_dir,"--cores":str(threads)}
            
            mergedOpts={**newOpts,**kwargs}
            #run trimgalore
            self.run_trimgalore(verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
            """
            running trim galore will create one file named <input>_trimmed.fq
            move these files to the specified out files
            """
            oldFile=os.path.join(out_dir,pu.get_file_basename(fq)+"_trimmed.fq")
            
            pe.move_file(oldFile,out_file)
            
            if not pu.check_files_exist(out_file):
                print("Trimgalore failed")
                return ("",)
            return (out_file,)
        
        
            
    def run_trimgalore(self,valid_args=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper for running trimgalore
        
        Parameters
        ----------
        valid_args: list
            List of valid args
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        
        kwargs: dict
            Options to pass to trimgalore (will override existing parameters)
            
        :return: Status of trimgalore command
        :rtype: bool
        """
      
        
        #create command to run
        trimgalore_cmd=['trim_galore']
        trimgalore_cmd.extend(pu.parse_unix_args(valid_args,kwargs))
        
        
        #start ececution
        status=pe.execute_command(trimgalore_cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            pu.print_boldred("trimgalore failed")
        
        #return status
        return status
        

            
            

class BBmap(RNASeqQC):
    """This class represents bbmap programs
    """
    def __init__(self,threads=None,max_memory=None):
        """
        Parameters
        ----------
        
        threads: int
            num threads to use
        max_memory: Max memory to use in GB
        """
        #run super to inherit parent class properties
        super().__init__() 
        self.programName="bbduk.sh"
        self.dep_list=[self.programName]
        #check if program exists
        if not pe.check_dependencies(self.dep_list):
            raise Exception("ERROR: "+ self.programName+" not found.")
            
        """
        self.valid_args=['in','in2','ref','literal','touppercase','interleaved','qin','reads','copyundefined',
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
        """
        
        #use max threads by default
        if not threads:
            threads=os.cpu_count()
        self.threads=threads
        
        #use floor(max available memory) by default
        if not max_memory:
            total_mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')  
            total_mem_gib = total_mem_bytes/(1024.**3)
            max_memory=math.floor(total_mem_gib)
        
        self.max_memory=max_memory
            
            
            
            
    def perform_qc(self,sra_object,out_dir="",out_suffix="_bbduk",overwrite=True,threads=None,max_memory=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Run bbduk on fastq files specified by the sra_object
        
        Parameters
        ----------
        
        sra_object: SRA
            an SRA object
        out_dir: string
            Path to out dir. Default: sra_object.location
        out_suffix: string
            Suffix for output file name
        overwrite: bool
            overwrite existing files
        threads: int
            Num threads to use
        max_memory: float
            Max memory to use in GB
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        
        kwargs: dict
            options passed to bbduk
            
        :return: Returns the path of fastq files after QC. tuple has one item for single end files and 2 for paired.
        :rtype: tuple
            
        """
        
        #make out_dir
        if not out_dir:
                out_dir=sra_object.location
        else:
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
                
        if not threads:
            threads=self.threads
        if not max_memory:
            max_memory=self.max_memory
            
        memory_flag="-Xmx"+str(max_memory)+"g"
        
        #optimize parameters
        #if optimize:
        #    print("generating suggested parameters XXX TD")
                    
        if sra_object.layout=='PAIRED':
            fq1=sra_object.localfastq1Path
            fq2=sra_object.localfastq2Path
            
            out_fileName1=pu.get_file_basename(fq1)+out_suffix+".fastq"
            out_fileName2=pu.get_file_basename(fq2)+out_suffix+".fastq"
            out_file1Path=os.path.join(out_dir,out_fileName1)
            out_file2Path=os.path.join(out_dir,out_fileName2)
            
            newOpts={"in":fq1,"in2":fq2,"out":out_file1Path,"out2":out_file2Path,"--":(memory_flag,),"threads":str(threads)}
            mergedOpts={**newOpts,**kwargs}
            
            #run bbduk
            if self.run_bbduk(verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts):
                if pu.check_files_exist(out_file1Path,out_file2Path):
                    return(out_file1Path,out_file2Path)
            return("",)
            
            
        else:
            fq=sra_object.localfastqPath
            out_fileName=pu.get_file_basename(fq)+out_suffix+".fastq"
            out_filePath=os.path.join(out_dir,out_fileName)
            newOpts={"in":fq,"out":out_filePath,"--":(memory_flag,)}
            mergedOpts={**newOpts,**kwargs}
            
            #run bbduk
            if self.run_bbduk(verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts):
                if pu.check_files_exist(out_filePath):
                    return(out_filePath,)
            return("",)
                    
    
    
    
    
    def run_bbduk(self,valid_args=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper to run bbduk.sh
        valid_args: list
            A list of valid arguments
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            options passed to bbduk
        
        """
              
        #create command to run
        bbduk_cmd=["bbduk.sh"]
        
        #bbduk.sh follows java style arguments
        bbduk_cmd.extend(pu.parse_java_args(valid_args,kwargs))
        
        
        #start ececution
        status=pe.execute_command(bbduk_cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            pu.print_boldred("bbduk failed")
        #return status
        return status
    
 
    
    def perform_cleaning(self,sra_object,bbsplit_index,out_dir="",out_suffix="_bbsplit",overwrite=True,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """
        Remove contaminated reads mapping to given reference using bbsplit
        
        Parameters
        ----------
        
        sra_object: SRA
            an SRA object
        bbsplit_index: string
            Path to bbsplit index or fasta file which will generate index
        out_dir: string
            Path to output dir. Default: sra_object.location
        out_suffix: string
            Suffix for output file name
        overwrite: bool
            overwrite existing files
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        
        kwargs: dict
            options passed to bbsplit

            :return: Returns the path of fastq files after QC. tuple has one item for single end files and 2 for paired.
            :rtype: tuple
        """
        
        #check index
        indexPath=""
        if not pu.check_paths_exist(bbsplit_index):
            #index folder doesn't exist
            #check if input is path to fasta
            if not pu.check_files_exist(bbsplit_index):
                print("Error: Please check bbsplit index")
                return ("",)
            #check if index folder "ref" exists in this directory
            indexPath=os.path.join(pu.get_file_directory(bbsplit_index),"ref")
            if pu.check_paths_exist(indexPath):
                print("Using bbsplit index: "+indexPath)
            else:
                #create new index
                print("Creating new index"+indexPath)
                newOpts={"ref_x":bbsplit_index,"path": pu.get_file_directory(bbsplit_index)}
                mergedOpts={**kwargs,**newOpts}
                #run bbduk
                if not self.run_bbsplit(verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts):
                    print("Error creating bbsplit index.")
                    return ("",)
                if not pu.check_paths_exist(indexPath):
                    print("Error creating bbsplit index.")
                    return ("",)
        else:
            indexPath=bbsplit_index
                
        
        #indexPath point to the ref directory, go one directory higher
        indexPath=os.path.dirname(indexPath)
        
        
        #make out_dir
        if not out_dir:
                out_dir=sra_object.location
        else:
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
        
        if sra_object.layout=='PAIRED':
            fq1=sra_object.localfastq1Path
            fq2=sra_object.localfastq2Path
            #append input and output options
            
            out_fileName1=pu.get_file_basename(fq1)+out_suffix+".fastq"
            out_fileName2=pu.get_file_basename(fq2)+out_suffix+".fastq"
            out_file1Path=os.path.join(out_dir,out_fileName1)
            out_file2Path=os.path.join(out_dir,out_fileName2)
            
            newOpts={"in1":fq1,"in2":fq2,"outu1":out_file1Path,"outu2":out_file2Path,"path":indexPath}
            mergedOpts={**kwargs,**newOpts}
            
            #run bbsplit
            if self.run_bbsplit(verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts):
                if pu.check_files_exist(out_file1Path,out_file2Path):
                    return(out_file1Path,out_file2Path)
            return("",)
            
            
        else:
            fq=sra_object.localfastqPath
            #append input and output options
           
            out_fileName=pu.get_file_basename(fq)+out_suffix+".fastq"
            out_filePath=os.path.join(out_dir,out_fileName)
            newOpts={"in":fq,"outu":out_filePath,"path":indexPath}
            mergedOpts={**kwargs,**newOpts}
            
            #run bbsplit
            if self.run_bbsplit(verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts):
                if pu.check_files_exist(out_filePath):
                    return(out_filePath,)
            
            return("",)
    
    
    
    def run_bbsplit(self,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """wrapper to run bbsplit
        
        :return: Status of bbsplit command
        :rtype: bool
        """
        
        bbsplit_args=['ref','ref_x','build','path','in','in1','in2','outu','outu2','outu1','qin','interleaved',
                          'maxindel','minratio','minhits','ambiguous','ambiguous2',
                          'qtrim','untrim','out_','basename','bs','scafstats',
                          'refstats','nzo','-Xmx','-eoom','-da']
        
        
        #create command to run
        bbsp_cmd=["bbsplit.sh"]
        
        #bbduk.sh follows java style arguments
        bbsp_cmd.extend(pu.parse_java_args(bbsplit_args,kwargs))
        
        
        #start ececution
        status=pe.execute_command(bbsp_cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            pu.print_boldred("bbsplit failed")
        #return status
        return status
    
    
    
    
    
    
        
