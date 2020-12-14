#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 19:53:42 2019

@author: usingh
contains classes of RNA-Seq alignment programs
"""
import os
from pyrpipe.runnable import Runnable
from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe
from pyrpipe import tools
from pyrpipe import valid_args
from pyrpipe import param_loader as pl
from pyrpipe import _dryrun
from pyrpipe import _threads
from pyrpipe import _params_dir


class Aligner(Runnable):
    """This is an abstract class for alignment programs.
    """
    def __init__(self,*args,index=None,genome=None,**kwargs):
        super().__init__(*args,**kwargs)
        self._category="Aligner"
        self._command=None
        self.index=index
        self.genome=genome
       
    def build_index(self):
        """function to create an index used by the aligner
        """
        pass    
    def check_index(self):
        """Function to check if index of this object is valid and exists
        """   
    def perform_alignment(self,sra_object):
        """Function to perform alignment taking and sraobject as input
        
        """
        pass

class Hisat2(Aligner):
    """
    Attributes
    ----------
    
    """ 
    def __init__(self,*args,index=None,genome=None,**kwargs):
        super().__init__(*args,**kwargs)
        self._command='hisat2'
        self._deps=[self._command,'samtools']
        self.index=index
        self.genome=genome
        self._param_yaml='hisat2.yaml'
        self._valid_args=valid_args._args_HISAT2
        self.check_dependency()
        self.init_parameters(*args,**kwargs)
        
        #check index flag in args
        if not index and '-x' in self._kwargs:
            self.index=self._kwargs['-x']
        
        #check index
        if not self.check_index():
            if not (pu.check_files_exist(self.genome)):
                pu.print_boldred("Hisat2 index '{}' not found; New index could not be created as genome file '{}' not found.".format(self.index,self.genome))
                raise Exception("Please provide a valid Hisat2 index, or a valid fasta file to generate the index")
            else:
                #call build index to generate index
                self.build_index(self.index,self.genome)
                #set index 
                self._kwargs['-x']=self.index
                

    def build_index(self,index_path,genome,overwrite=False,objectid="NA"):
        """Build a hisat index with given parameters and saves the new index to self.index.
        
        Parameters
        ----------
        
        index_path: string
            Path where the index will be created
            
        index_name: string
            A name for the index
            
        args: tuple
            Path to reference input files
        threads: int
            Num threads to use
            
        verbose : bool
            Print stdout and std error
            
        quiet : bool 
            Print nothing
            
        logs : bool 
            Log this command to pyrpipe logs
            
        objectid : string 
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        
        kwargs: dict
            Parameters for the hisat2-build command
        
        
            
        :return: Returns the status of hisat2-build
        :rtype: bool
        """
        
        #if index already exists then exit
        if not overwrite:
            #check if files exists
            if pu.check_hisatindex(index_path):
                pu.print_green("Hisat2 index {} already exists.".format(index_path))
                self.index=os.path.join(index_path)
                return True
            
        #check input files
        if not pu.check_files_exist(genome):
            pu.print_boldred("Please provide a valid input fasta file to build Hisat2 index")
            raise Exception("Please check input to hisat2 build index")
            return False
        
        indexdir=pu.get_file_directory(index_path)
        #create the out dir
        if not pu.check_paths_exist(indexdir):
            if not pu.mkdir(indexdir):
                raise Exception("Error creating hisat2 index. Failed to create index directory.")
                
        
        hisat2Buildvalid_args=['-c','--large-index','-a','-p','--bmax','--bmaxdivn','--dcv','--nodc','-r','-3','-o',
                                  '-t','--localoffrate','--localftabchars','--snp','--haplotype','--ss','--exon',
                                  '--seed','-q','-h','--usage','--version']
        
        args=(genome,index_path)
        internal_kwargs={"-p":_threads}
        #read build parameters
        yamlfile=os.path.join(_params_dir,'hisat2_index.yaml')
        if pu.check_files_exist(yamlfile):
            yaml_params=pl.YAML_loader(yamlfile)
            yaml_kwargs=yaml_params.get_kwargs()
            internal_kwargs={**yaml_kwargs,**internal_kwargs}
            
        #add positional args
        internal_kwargs['--']=args
        
        hisat2Build_Cmd=['hisat2-build']
        hisat2Build_Cmd.extend(pu.parse_unix_args(hisat2Buildvalid_args,internal_kwargs))
        
        #execute command
        status=pe.execute_command(hisat2Build_Cmd,objectid=objectid)
        
        if status:
            if pu.check_hisatindex(index_path) and not _dryrun:
                #update object's index
                self.index=index_path
                if self.check_index():
                    return True
        else:
            raise Exception("Error building Hisat2 index")
        
        return True
        
    def perform_alignment(self,sra_object,out_suffix="_hisat2",out_dir="",objectid="NA"):
        """Function to perform alignment using sra_object.
        
        Parameters
        ----------
        
        sra_object SRA object
            An object of type SRA. The path to fastq files will be obtained from this object.
        out_suffix: string
            Suffix for the output sam file
        threads: int
            Num threads to use
        overwrite: bool
            Overwrite output sam if already exist
        verbose: bool
            Print stdout and std error
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to hisat2. This will override the existing options in self.passed_args_dict (only replace existing arguments and not replace all the arguments).
        """
        #check out dir
        if not out_dir:
            out_dir=sra_object.directory
        else:
            #create out_dir if not exists
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
            
        #create path to output sam file
        outSamFile=os.path.join(out_dir,sra_object.srr_accession+out_suffix+".sam")
        outBamFile=os.path.join(out_dir,sra_object.srr_accession+out_suffix+"_sorted.bam")
        
 
        #find layout and fq file paths
        if sra_object.layout == 'PAIRED':
            internal_kwargs={"-1":sra_object.fastq_path,"-2":sra_object.fastq2_path,"-S":outSamFile,"-p":_threads,"-x":self.index}
        else:
            internal_kwargs={"-U":sra_object.fastq_path,"-S":outSamFile,"-p":_threads,"-x":self.index}
        
        #call run_hisat2
        status=self.run(None,objectid=sra_object.srr_accession,target=outBamFile,**internal_kwargs)
        
        if status:
            if not pu.check_files_exist(outSamFile) and not _dryrun:
                return ""
            #convert to bam before returning
            return tools.Samtools().sam_sorted_bam(outSamFile)
            #return outSamFile
        
        return ""
            
     
    
    def check_index(self):
        if hasattr(self,'index'):
            return(pu.check_hisatindex(self.index))
        else:
            return False


class Star(Aligner):
    """This class represents STAR program.
    
    Attributes
    ----------
    """ 
    def __init__(self,*args,index=None,genome=None,**kwargs):
        super().__init__(*args,**kwargs)
        self._command='STAR'
        self._deps=[self._command]
        self.index=index
        self.genome=genome
        self._param_yaml='star.yaml'
        self._valid_args=valid_args._args_STAR
        self.check_dependency()
        self.init_parameters(*args,**kwargs)
               
        #check index flag in args
        if not index and '--genomeDir' in self._kwargs:
            self.index=self._kwargs['--genomeDir']
            
        #check index
        if not pu.check_starindex(index):
            if not (pu.check_files_exist(self.genome)):
                pu.print_boldred("STAR index '{}' not found; New index could not be created as genome file '{}' not found.".format(self.index,self.genome))
                raise Exception("Please provide a valid STAR index, or a valid fasta file to generate the index")
            else:
                #call build index to generate index
                self.build_index(index,self.genome)
                #set index 
                self._kwargs['--genomeDir']=self.index
    
    def build_index(self,index_path,genome,overwrite=False,objectid="NA"):
        """Build a star index with given parameters and saves the new index to self.index.
        
        Parameters
        ----------
        
        index_path: string
            Path where the index will be created
        genome: tuple
            Path to reference input files
        overwrite: bool
            Overwrite if index already exists
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
            
        :return: Returns status of star command
        :rtype: bool
        """
        
        #if index already exists then exit
        if not overwrite:
            if pu.check_starindex(index_path):
                pu.print_green("STAR index {} already exists.".format(index_path))
                self.index=index_path
                return True
            
        #check input files
        if not (pu.check_files_exist(genome)):
            pu.print_boldred("Please provide a valid input fasta file to build STAR index")
            raise Exception("Please check input to build star index")
            
        
        #create index path if doesnt exist
        if not pu.check_paths_exist(index_path):
            if not pu.mkdir(index_path):
                raise Exception("Error creating STAR index. Failed to create index directory.")
                return False
        
        
        #determine parameters and execute cmd
        #internal_args=()
        internal_kwargs={"--runMode":"genomeGenerate","--genomeDir":index_path,"--genomeFastaFiles":genome,"--runThreadN":_threads}
        
        #read build parameters
        yamlfile=os.path.join(_params_dir,'star_index.yaml')
        if pu.check_files_exist(yamlfile):
            yaml_params=pl.YAML_loader(yamlfile)
            yaml_kwargs=yaml_params.get_kwargs()
            internal_kwargs={**yaml_kwargs,**internal_kwargs}
        
        starbuild_Cmd=['STAR']
        starbuild_Cmd.extend(pu.parse_unix_args(valid_args._args_STAR,internal_kwargs))
        
        #execute command
        status=pe.execute_command(starbuild_Cmd,objectid=objectid)
        if status:
            if pu.check_paths_exist(index_path) and not _dryrun:
                #update object's index
                self.index=index_path
                if self.check_index():
                    return True
        else:
            raise Exception("Error building STAR index")
        
        return True
            
    def perform_alignment(self,sra_object,out_suffix="_star",out_dir="",objectid="NA"):
        """Function to perform alignment using star and the provided SRA object.
        All star output will be written to the sra_object directory by default.
        
        Parameters
        ----------
        
        sra_object: SRA object
            An object of type SRA. The path to fastq files will be obtained from this object.
        out_suffix: string
            Suffix for the output file
        out_dir: str
            outout directory default: sra_object.directory
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
        out_type: str
            Out type options for star: sam, sorted_bam, unsorted_bam [Default: sorted_bam]
        optimize: bool
            If true optimize maping parameters based on read length [Default: True]
        threads: int
            Num threads to use. If supplied will override threads supplied during __init__
        kwargs: dict
            Options to pass to STAR. This will override the existing options in self.passed_args_dict (only replace existing arguments and not replace all the arguments).
        

        :return: Return the path to the output dir
        :rtype: string
        """        
        
        if not out_dir:
            out_dir=sra_object.directory
        else:
            #create out_dir if not exists
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
        
        #find layout and fq file paths
        if sra_object.layout == 'PAIRED':
            internal_kwargs={"--readFilesIn":sra_object.fastq_path+" "+sra_object.fastq2_path}
        else:
            internal_kwargs={"--readFilesIn":sra_object.fastq_path}
        #add out dir
        internal_kwargs["--outFileNamePrefix"]=out_dir+"/"
        #threads
        internal_kwargs["--runThreadN"]=_threads
        #add index
        internal_kwargs["--genomeDir"]=self.index
        
        #the expected out file
        bam=os.path.join(out_dir,'Aligned.out.bam')
        #call star
        status=self.run(None,objectid=sra_object.srr_accession,target=bam,**internal_kwargs)
                
        
        if status:
            #return rename the bam  file and return path
            #star can return Aligned.sortedByCoord.out.bam Aligned.out.bam Aligned.toTranscriptome.out.bam  
            #return sorted bam or unsorted bam which ever is present
            if 'SortedByCoordinate' in self._kwargs['--outSAMtype']:
                bam=os.path.join(out_dir,'Aligned.sortedByCoord.out.bam')
                
            finalbam=bam.split('.bam')[0]+out_suffix+'.bam'
            
            if not _dryrun:
                pe.move_file(bam,finalbam)
                if not pu.check_files_exist(finalbam):
                    return ""
                
            return finalbam
        
        return ""
 
    def check_index(self):
        if hasattr(self,'index'):
            return(pu.check_starindex(self.index))
        else:
            return False
            

class Bowtie2(Aligner):
    """This class represents bowtie2 program.
    
       Parameters
       ----------
       
       index: string
            path to a bowtie2 index. This index will be used when bowtie2 is invoked using this object.
       threads: int
            Num threads to use
            
       Attributes
       ----------
    """ 
    def __init__(self,*args,index=None,genome=None,**kwargs):
        """Bowtie2 constructor. Initialize bowtie2 index and other parameters.
        """       
        super().__init__(*args,**kwargs)
        self._command='bowtie2'
        self._deps=[self._command]
        self.index=index
        self.genome=genome
        self._param_yaml='bowtie2.yaml'
        self._valid_args=valid_args._args_BOWTIE2
        self.check_dependency()
        self.init_parameters(*args,**kwargs)
               
        #check index flag in args
        if not index and '-x' in self._kwargs:
            self.index=self._kwargs['-x']
        
        #if index is passed, update the passed arguments
        if not pu.check_bowtie2index(index):
            if not (pu.check_files_exist(self.genome)):
                pu.print_boldred("Bowtie2 index '{}' not found; New index could not be created as genome file '{}' not found.".format(self.index,self.genome))
                raise Exception("Please provide a valid Bowtie2 index, or a valid fasta file to generate the index")
            else:
                #call build index to generate index
                self.build_index(index,self.genome)
                #set index 
                self._kwargs['-x']=self.index
            
        
    def build_index(self,index_path,genome,overwrite=False,objectid="NA"):
        """Build a bowtie2 index with given parameters and saves the new index to self.index.
        
        Parameters
        ----------
        
        index_path: string
            Path where the index will be created
        index_name: string
            A name for the index
        args: tuple
            Path to reference input files
        threads: int
            Num threads to use
        overwrite: bool
            Overwrite already existing index
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to bowtie2. This will override the existing options in self.passed_args_dict (only replace existing arguments and not replace all the arguments).
            
        
 
        :return: Returns the status of bowtie-build
        :rtype: bool
        """
        
        #check input references
        if not overwrite:
            if pu.check_bowtie2index(index_path):
                pu.print_green("bowtie index {} already exists.".format(index_path))
                self.index=index_path
                return True
        
        #check input files
        if not (pu.check_files_exist(genome)):
            pu.print_boldred("Please provide a valid input fasta file to build bowtie2 index")
            raise Exception("Please check input to star build index")
            return False
        
        
        bowtie2_build_args=['-f','-c','--large-index','--debug','--sanitized','--verbose','-a',
                            '--noauto','-p','--packed','--bmax','--bmaxdivn','--dcv','--nodc',
                            '-r','--noref','-3','--justref','-o','--offrate','-t','--ftabchars',
                            '--threads','--seed','-q','--quiet']
        
        #create the out dir
        indexdir=pu.get_file_directory(index_path)
        if not pu.check_paths_exist(indexdir):
            if not pu.mkdir(indexdir):
                raise Exception("Error creating bowtie2 index. Failed to create index directory.")
                return False
        
        args=(genome,index_path)
        internal_kwargs={"--threads":_threads}
        
        #read build parameters
        yamlfile=os.path.join(_params_dir,'bowtie2_index.yaml')
        if pu.check_files_exist(yamlfile):
            yaml_params=pl.YAML_loader(yamlfile)
            yaml_kwargs=yaml_params.get_kwargs()
            internal_kwargs={**yaml_kwargs,**internal_kwargs}
        
        #add positional args
        internal_kwargs['--']=args
        
        bowtie2Build_Cmd=['bowtie2-build']       
        #add options
        bowtie2Build_Cmd.extend(pu.parse_unix_args(bowtie2_build_args,internal_kwargs))
    
        
        #start ececution
        status=pe.execute_command(bowtie2Build_Cmd,objectid=objectid)
        if not status:
            pu.print_boldred("bowtie2-build failed")
            return False
        
        if status:
            if pu.check_bowtie2index(index_path) and not _dryrun:
                #update object's index
                self.index=index_path
                if self.check_index():
                    return True
        else:
            raise Exception("Error building bowtie2 index")
            
        
        return True
        
        
    
    def perform_alignment(self,sra_object,out_suffix="_bowtie2",out_dir="",objectid="NA"):
        """Function to perform alignment using self object and the provided sra_object.
        
        Parameters
        ----------
        
        sra_object: SRA object
            An object of type SRA. The path to fastq files will be obtained from this object.
        out_suffix: string
            Suffix for the output sam file
        out_dir: str
            Path to out dir
        threads: int
            Num threads to use
        overwrite: bool
            Overwrite sam file if exixts
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to bowtie2. This will override the existing options in self.passed_args_dict (only replace existing arguments and not replace all the arguments).
        
        :return: Returns the output sam file path
        :rtype: string
        """
        if not out_dir:
            out_dir=sra_object.directory
        else:
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
                
        #create path to output sam file
        outSamFile=os.path.join(out_dir,sra_object.srr_accession+out_suffix+".sam")
        outBamFile=os.path.join(out_dir,sra_object.srr_accession+out_suffix+"_sorted.bam")
                    
        
        #find layout and fq file paths
        if sra_object.layout == 'PAIRED':
            internal_kwargs={"-1":sra_object.fastq_path,"-2":sra_object.fastq2_path,"-S":outSamFile,"--threads":_threads,"-x":self.index}
        else:
            internal_kwargs={"-U":sra_object.fastq_path,"-S":outSamFile,"--threads":_threads,"-x":self.index}
        
        
        
        status=self.run(None,objectid=sra_object.srr_accession,target=outBamFile,**internal_kwargs)
        
        
        if status:
            if not pu.check_files_exist(outSamFile) and not _dryrun:
                return ""
            #convert to bam before returning
            return tools.Samtools().sam_sorted_bam(outSamFile)
            #return outSamFile
        
        return ""
        
        return ""
    
    def check_index(self):
        """Function to check bowtie index.
        Returns True is index exist on disk.
        """
        if hasattr(self,'index'):
            return(pu.check_bowtie2index(self.index))
        return False



         
