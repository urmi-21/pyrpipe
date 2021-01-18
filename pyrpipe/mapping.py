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
from pyrpipe import _force
from pyrpipe import _params_dir



class Aligner(Runnable):
    """This is an abstract class for alignment programs.
    """
    def __init__(self,*args,index=None,genome=None,threads=None,**kwargs):
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
    Extends Aligner class
    Attributes
    ----------
    
    """ 
    def __init__(self,*args,index=None,genome=None,threads=None,**kwargs):
        """
        init Hisat2 object

        Parameters
        ----------
        *args : tuple
            Positional arguements
        index : Str, optional
            Path to Hisat index. If index is not present it will generate an index using the genome. Index can be supplied via the hisat2.yaml file too. The default is None.
        genome : Str, optional
            Path to the reference genome. This will be used to build an index if index is not present The default is None.
        threads : int, optional
            Threads to use for hisat2. This will override the global --threads parameter supplied to pyrpipe. The default is None.
        **kwargs : dict
            keyword arguments

        Raises
        ------
        ValueError
            Raises ValueError if hisat index is not found and genome is not present to generate an index.

        Returns
        -------
        None.

        """
        super().__init__(*args,**kwargs)
        self._command='hisat2'
        self._deps=[self._command,'samtools']
        self.index=index
        self.genome=genome
        self._param_yaml='hisat2.yaml'
        self._valid_args=valid_args._args_HISAT2       
        
        #resolve threads to use
        self.resolve_parameter("-p",threads,_threads,'_threads')
        #resolve index to use
        self.resolve_parameter("-x",index,index,'index')
        
        #check index
        if not self.check_index():
            if not (pu.check_files_exist(self.genome)):
                pu.print_boldred("Hisat2 index '{}' not found; New index could not be created as genome file '{}' not found.".format(self.index,self.genome))
                raise ValueError("Please provide a valid Hisat2 index, or a valid fasta file to generate the index")
            else:
                #call build index to generate index
                self.build_index(self.index,self.genome)
                #set index 
                #self._kwargs['-x']=self.index
        
                

    def build_index(self,index_path,genome,objectid="NA"):
        """Build a hisat index with given parameters and saves the new index to self.index.
        
        Parameters
        ----------
        
        index_path: string
            Path where the index will be created
        genome: string
            Path to the reference genome
        objectid : string 
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
            
        :return: Returns the status of hisat2-build
        :rtype: bool
        """
        
        #if index already exists then exit
        if not _force:
            #check if files exists
            if pu.check_hisatindex(index_path):
                pu.print_green("Hisat2 index {} already exists.".format(index_path))
                self.index=os.path.join(index_path)
                return True
            
        #check input files
        if not pu.check_files_exist(genome):
            pu.print_boldred("Please provide a valid input fasta file to build Hisat2 index")
            raise ValueError("Please check input to hisat2 build index")
        
        indexdir=pu.get_file_directory(index_path)
        #create the out dir
        if not pu.check_paths_exist(indexdir):
            if not pu.mkdir(indexdir):
                raise OSError("Error creating hisat2 index. Failed to create index directory.")
                
        
        hisat2Buildvalid_args=valid_args._args_HISAT2BUILD
        
        
        args=(genome,index_path)
        internal_kwargs={"-p":self._threads}
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
            raise OSError("Error building Hisat2 index")
        
        return True
        
    def perform_alignment(self,sra_object,out_suffix="_hisat2",out_dir="",objectid="NA"):
        """Function to perform alignment using sra_object.
        
        Parameters
        ----------
        
        sra_object SRA object
            An object of type SRA. The path to fastq files will be obtained from this object.
        out_suffix: string
            Suffix for the output sam file
        out_dir: string
            Directory to save the results. Default value is sra_object.directory
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        
        :return: Returns the sorted bam file path after converting sam to bam and sorting it
        :rtype: string
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
        #check if final bam already exists
        if not _force and pu.check_files_exist(outBamFile):
            pu.print_green('Target files {} already exist.'.format(outBamFile))
            return outBamFile
        
 
        #find layout and fq file paths
        if sra_object.layout == 'PAIRED':
            internal_kwargs={"-1":sra_object.fastq_path,"-2":sra_object.fastq2_path,"-S":outSamFile}
        else:
            internal_kwargs={"-U":sra_object.fastq_path,"-S":outSamFile}
        
        #call run_hisat2
        status=self.run(None,objectid=sra_object.srr_accession,target=outSamFile,**internal_kwargs)
        
        if status:
            if not pu.check_files_exist(outSamFile) and not _dryrun:
                return ""
            #convert to bam before returning; returns outBamFile
            return tools.Samtools().sam_sorted_bam(outSamFile)
            #return outSamFile
        
        return ""
     
    
    def check_index(self):
        """
        Check self.index exists

        Returns
        -------
        bool
            True if index exists.

        """
        if hasattr(self,'index'):
            return(pu.check_hisatindex(self.index))
        else:
            return False


class Star(Aligner):
    """This class represents STAR program.
    Extends the Aligner class
    Attributes
    ----------
    """ 
    def __init__(self,*args,index=None,genome=None,threads=None,**kwargs):
        """
        init Star object

        Parameters
        ----------
        *args : tuple
            Positional arguements
        index : Str, optional
            Path to Star index. If index is not present it will generate an index using the genome. Index can be supplied via the star.yaml file too. The default is None.
        genome : Str, optional
            Path to the reference genome. This will be used to build an index if index is not present. The default is None.
        threads : int, optional
            Threads to use for STAR. This will override the global --threads parameter supplied to pyrpipe. The default is None.
        **kwargs : dict
            keyword arguments

        Raises
        ------
        ValueError
            Raises ValueError if STAR index is not found and genome is not present to generate an index.

        Returns
        -------
        None.

        """
        super().__init__(*args,**kwargs)
        self._command='STAR'
        self._deps=[self._command]
        self.index=index
        self.genome=genome
        self._param_yaml='star.yaml'
        self._valid_args=valid_args._args_STAR
        
        #resolve threads to use
        self.resolve_parameter("--runThreadN",threads,_threads,'_threads')
        #resolve index to use
        self.resolve_parameter("--genomeDir",index,index,'index')
            
        #check index
        if not self.check_index():
            if not (pu.check_files_exist(self.genome)):
                pu.print_boldred("STAR index '{}' not found; New index could not be created as genome file '{}' not found.".format(self.index,self.genome))
                raise ValueError("Please provide a valid STAR index, or a valid fasta file to generate the index")
            else:
                #call build index to generate index
                self.build_index(index,self.genome)

    
    def build_index(self,index_path,genome,objectid="NA"):
        """Build a STAR index with given parameters and saves the new index to self.index.
        
        Parameters
        ----------
        
        index_path: string
            Path where the index will be created
        genome: string
            Path to the reference genome
        objectid : string 
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
            
        :return: Returns the status of STAR-build index
        :rtype: bool
        """
        
        
        #if index already exists then exit
        if not _force:
            if pu.check_starindex(index_path):
                pu.print_green("STAR index {} already exists.".format(index_path))
                self.index=index_path
                return True
            
        #check input files
        if not (pu.check_files_exist(genome)):
            pu.print_boldred("Please provide a valid input fasta file to build STAR index")
            raise ValueError("Please check input to build star index")
            
        
        #create index path if doesnt exist
        if not pu.check_paths_exist(index_path):
            if not pu.mkdir(index_path):
                raise OSError("Error creating STAR index. Failed to create index directory.")
                return False
        
        
        #determine parameters and execute cmd
        #internal_args=()
        internal_kwargs={"--runMode":"genomeGenerate","--genomeDir":index_path,"--genomeFastaFiles":genome,"--runThreadN":self._threads}
        
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
            raise OSError("Error building STAR index")
        
        return True
            
    def perform_alignment(self,sra_object,out_suffix="_star",out_dir="",objectid="NA"):
        """Function to perform STAR alignment using sra_object.
        
        Parameters
        ----------
        
        sra_object SRA object
            An object of type SRA. The path to fastq files will be obtained from this object.
        out_suffix: string
            Suffix for the output sam file
        out_dir: string
            Directory to save the results. Default value is sra_object.directory
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        
        :return: Returns the path to output bam
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
       
        
        #the expected out file
        #star can return Aligned.sortedByCoord.out.bam Aligned.out.bam Aligned.toTranscriptome.out.bam  
        #return sorted bam or unsorted bam which ever is present
        bam=os.path.join(out_dir,'Aligned.out.bam')
        
        #if outSAMtype is not specified make it bam by default
        if not '--outSAMtype' in self._kwargs: self._kwargs['--outSAMtype']='BAM SortedByCoordinate'
        
        if '--outSAMtype' in self._kwargs and 'SortedByCoordinate' in self._kwargs['--outSAMtype']:
                bam=os.path.join(out_dir,'Aligned.sortedByCoord.out.bam')
        finalbam=bam.split('.bam')[0]+out_suffix+'.bam'
        
        #check if final bam already exists
        if not _force and pu.check_files_exist(finalbam):
            pu.print_green('Target files {} already exist.'.format(finalbam))
            return finalbam
            
        #call star
        status=self.run(None,objectid=sra_object.srr_accession,target=bam,**internal_kwargs)
                
        
        if status:
            #return rename the bam  file and return path
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
    """This Bowtie2 aligner class.
       Extends Aligner class
       
       Attributes
       ----------
    """ 
    def __init__(self,*args,index=None,genome=None,threads=None,**kwargs):
        """
        init Bowtie2 object

        Parameters
        ----------
        *args : tuple
            Positional arguements
        index : Str, optional
            Path to Bowtie2 index. If index is not present it will generate an index using the genome. Index can be supplied via the bowtie2.yaml file too. The default is None.
        genome : Str, optional
            Path to the reference genome. This will be used to build an index if index is not present The default is None.
        threads : int, optional
            Threads to use for Bowtie2. This will override the global --threads parameter supplied to pyrpipe. The default is None.
        **kwargs : dict
            keyword arguments

        Raises
        ------
        ValueError
            Raises ValueError if hista index is not found and genome is not present to generate an index.

        Returns
        -------
        None.

        """     
        super().__init__(*args,**kwargs)
        self._command='bowtie2'
        self._deps=[self._command]
        self.index=index
        self.genome=genome
        self._param_yaml='bowtie2.yaml'
        self._valid_args=valid_args._args_BOWTIE2
        
        
        #resolve threads to use
        self.resolve_parameter("-p",threads,_threads,'_threads')
        #resolve index to use
        self.resolve_parameter("-x",index,index,'index')
        
        #if index is passed, update the passed arguments
        if not self.check_index():
            if not (pu.check_files_exist(self.genome)):
                pu.print_boldred("Bowtie2 index '{}' not found; New index could not be created as genome file '{}' not found.".format(self.index,self.genome))
                raise ValueError("Please provide a valid Bowtie2 index, or a valid fasta file to generate the index")
            else:
                #call build index to generate index
                self.build_index(self.index,self.genome)
            
        
    def build_index(self,index_path,genome,objectid="NA"):
        """Build a bowtie2 index with given parameters and saves the new index to self.index.
        
        Parameters
        ----------
        
        index_path: string
            Path where the index will be created
        genome: string
            Path to the reference genome
        objectid : string 
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
            
        :return: Returns the status of bowtie2-build
        :rtype: bool
        """
        
        #check input references
        if not _force:
            if pu.check_bowtie2index(index_path):
                pu.print_green("bowtie index {} already exists.".format(index_path))
                self.index=index_path
                return True
        
        #check input files
        if not (pu.check_files_exist(genome)):
            pu.print_boldred("Please provide a valid input fasta file to build bowtie2 index")
            raise ValueError("Please check input to star build index")
            return False
        
        
        bowtie2_build_args=['-f','-c','--large-index','--debug','--sanitized','--verbose','-a',
                            '--noauto','-p','--packed','--bmax','--bmaxdivn','--dcv','--nodc',
                            '-r','--noref','-3','--justref','-o','--offrate','-t','--ftabchars',
                            '--threads','--seed','-q','--quiet']
        
        #create the out dir
        indexdir=pu.get_file_directory(index_path)
        if not pu.check_paths_exist(indexdir):
            if not pu.mkdir(indexdir):
                raise OSError("Error creating bowtie2 index. Failed to create index directory.")
                return False
        
        args=(genome,index_path)
        internal_kwargs={"--threads":self._threads}
        
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
            raise OSError("Error building bowtie2 index")
            
        
        return True
        
        
    
    def perform_alignment(self,sra_object,out_suffix="_bowtie2",out_dir="",objectid="NA"):
        """Function to perform alignment using sra_object.
        
        Parameters
        ----------
        
        sra_object SRA object
            An object of type SRA. The path to fastq files will be obtained from this object.
        out_suffix: string
            Suffix for the output sam file
        out_dir: string
            Directory to save the results. Default value is sra_object.directory
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        :return: Returns the sorted bam file path after converting sam to bam and sorting it
        :rtype: string
        """
        if not out_dir:
            out_dir=sra_object.directory
        else:
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
                
        #create path to output sam file
        outSamFile=os.path.join(out_dir,sra_object.srr_accession+out_suffix+".sam")
        #outBamFile=os.path.join(out_dir,sra_object.srr_accession+out_suffix+"_sorted.bam")
                    
        
        #find layout and fq file paths
        if sra_object.layout == 'PAIRED':
            internal_kwargs={"-1":sra_object.fastq_path,"-2":sra_object.fastq2_path,"-S":outSamFile}
        else:
            internal_kwargs={"-U":sra_object.fastq_path,"-S":outSamFile}
        
        
        
        status=self.run(None,objectid=sra_object.srr_accession,target=outSamFile,**internal_kwargs)
        
        
        if status:
            if not pu.check_files_exist(outSamFile) and not _dryrun:
                return ""
            #convert to bam before returning; returns outBamFile
            return tools.Samtools().sam_sorted_bam(outSamFile)
            
        
        return ""
    
    
    
    def check_index(self):
        """Function to check bowtie index.
        Returns True is index exist on disk.
        """
        if hasattr(self,'index'):
            return(pu.check_bowtie2index(self.index))
        return False



         
