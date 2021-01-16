#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:21:01 2019

@author: usingh
"""
from pyrpipe.runnable import Runnable
from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe
import os
from pyrpipe import valid_args
from pyrpipe import _dryrun
from pyrpipe import _threads
from pyrpipe import _force


class Assembly(Runnable):
    """This class represents an abstract parent class for all programs which can perfrom transcripts assembly.
    """
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self._category="Assembler"
        self._command=None
        
    def perform_assembly(bam_file):
        """Function to perform assembly using a bam file as input. Inherited by all children.
        :param bam_file: path to input BAM
        :type bam_file: string
        
        :return: path to output GTF or output directory depending on the specific assembly program.
        :rtype: string
        """
        pass

class Stringtie(Assembly):
    """This class represents Stringtie program for transcript assembly.
        
        Parameters
        ----------
        
        threads: int
            number of threads
        guide: str
            Reference annotation gtf/gff to use as guide
            
        """
    def __init__(self,*args,threads=None,guide=None,**kwargs):
        super().__init__(*args,**kwargs)
        self._command='stringtie'
        self._deps=[self._command]
        self._param_yaml='stringtie.yaml'
        self._valid_args=valid_args._args_STRINGTIE
        
        #resolve threads to use
        self.resolve_parameter("-p",threads,_threads,'_threads')
        self.resolve_parameter("-G",guide,None,'_guide')
        

                         
    def perform_assembly(self,bam_file,out_dir=None,out_suffix="_stringtie",objectid="NA"):
        """Function to run stringtie using a bam file.
                
        Parameters
        ----------
        
        bam_file: string
            path to the bam file
        out_dir: string
            Path to out file
        out_suffix: string
            Suffix for the output gtf file
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        :return: Returns the path to output GTF file
        :rtype: string
        """
        
        #create path to output file
        fname=pu.get_file_basename(bam_file)
        
        if not out_dir:
            out_dir=pu.get_file_directory(bam_file)
        
        if not pu.check_paths_exist(out_dir):
            pu.mkdir(out_dir)
            
        out_gtf_file=os.path.join(out_dir,fname+out_suffix+".gtf")

        #Add output file name and input bam
        internal_args=(bam_file,)
        internal_kwargs={"-o":out_gtf_file}
        #add positional args
        internal_kwargs['--']=internal_args

        #call stringtie
        status=self.run(None,objectid=objectid,target=out_gtf_file,**internal_kwargs)
        
        if status:
            #check if sam file is present in the location directory of sraOb
            if not pu.check_files_exist(out_gtf_file) and not _dryrun:
                return ""
            return out_gtf_file
        
        return ""
    
    
class Cufflinks(Assembly):
    """This class represents cufflinks
    
    threads: int
            Number of threads to use
    """
    def __init__(self,*args,threads=None,guide=None,**kwargs):
        super().__init__(*args,**kwargs)
        self._command='cufflinks'
        self._deps=[self._command]
        self._param_yaml='cufflinks.yaml'
        self._valid_args=valid_args._args_CUFFLINKS
        
            
        #resolve threads to use
        self.resolve_parameter("-p",threads,_threads,'_threads')
        self.resolve_parameter("-G",guide,None,'_guide')
    
    def perform_assembly(self,bam_file,out_dir=None,out_suffix="_cufflinks",objectid="NA"):
        """Function to run cufflinks with BAM file as input.
                
        Parameters
        ----------
        bam_file: string
            path to bam file
        out_dir: 
            output directory
        out_suffix: string
            Suffix for the output gtf file
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession.
            
        :return: Returns the path to output GTF file
        :rtype: string       
        """
        
        #create path to output file
        fname=pu.get_file_basename(bam_file)
        if not out_dir:
            out_dir=pu.get_file_directory(bam_file)
        else:
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
        
            
        #Add output file name and input bam
        internal_args=(bam_file,)
        internal_kwargs={"-o":out_dir}
        #add positional args
        internal_kwargs['--']=internal_args
        
        #targets
        outfile=os.path.join(out_dir,"transcripts.gtf")
        out_gtf_file=os.path.join(out_dir,fname+out_suffix+".gtf")
        
        #if final file already exists
        if not _force and pu.check_files_exist(out_gtf_file):
            pu.print_green('Target files {} already exist.'.format(out_gtf_file))
            return out_gtf_file
        
        
        #call cufflinks
        status=self.run(None,objectid=objectid,target=outfile,**internal_kwargs)
        
        
        if status:
            if not _dryrun:
                pe.move_file(outfile,out_gtf_file)
                if not pu.check_files_exist(out_gtf_file):
                    return ""
                
            return out_gtf_file
            
            
        return ""
    
