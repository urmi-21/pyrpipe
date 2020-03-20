#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 18:20:41 2020

@author: usingh
"""

from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe
import os

class Quant:
    """This is an abstract class for quantification programs.
    """
    def __init__(self,index=""):
        self.category="Quantification"
        self.passedArgumentDict={}
        self.index=index
        
    def build_index(self):
        """function to create an index used by the quantification program
        """
        pass
    
    def check_index(self):
        """Function to check if index of this object is valid and exists
        """
    
    def perform_quant(self,sra_object):
        """Function to perform quant taking and sraobject as input
        
        """
        pass
    
class Kallisto(Quant):
    """This class represents kallisto
    
    kallisto_index: string
        path to kallisto index
    threads: int
        num threads to use
    """
    
    def __init__(self,kallisto_index,threads=None):
        super().__init__() 
        self.programName="kallisto"
        self.dep_list=[self.programName]        
        if not pe.check_dependencies(self.dep_list):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        """
        ##kallisto index
        self.validArgsIndex=['-i','--index','-k','--kmer-size','--make-unique']
        ##kallisto quant
        self.validArgsQuant=['-i','--index','-o','--output-dir','--bias','-b','--bootstrap-samples',
                             '--seed','--plaintext','--fusion','--single','--fr-stranded','--rf-stranded',
                             '-l','--fragment-length','-s','--sd','-t','--threads','--pseudobam']
        ##kallisto pseudo
        self.validArgsPseudo=['-i','--index','-o','--output-dir','-u','--umi','-b','--batch',
                              '--single','-l','--fragment-length','-s','--sd','-t','--threads','--pseudobam']
            ##kallisto h5dump
        self.validArgsh5dump=['-o','--output-dir']
        
        self.valid_args=pu.get_union(self.validArgsIndex,self.validArgsQuant,self.validArgsPseudo,self.validArgsh5dump)
        """
        
        if not threads:
            threads=os.cpu_count()
        self.threads=threads
        
        #if index is passed, update the passed arguments
        if len(kallisto_index)>0 and pu.check_files_exist(kallisto_index):
            print("kallisto index is: "+kallisto_index)
            self.kallisto_index=kallisto_index
        else:
            print("No kallisto index provided. Please use build_index() now to generate an index...")
            
    def build_index(self,index_path,index_name,fasta,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Function to  build kallisto index
        
        index_path: str
            path to the output directory
        index_name: str
            index name
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to kallisto. This will override the existing options in self.passed_args_dict (only replace existing arguments and not replace all the arguments).
            
        :return: Status of kallisto index
        :rtype: bool
        """
        
        #check input
        if not pu.check_files_exist(fasta):
            pu.print_boldred("{} does not exist. Exiting".format(fasta))
            return False
        
        #create out dir
        if not pu.check_paths_exist(index_path):
            if not pu.mkdir(index_path):
                print("ERROR in building kallisto index. Failed to create index directory.")
                return False
            
        indexOut=os.path.join(index_path,index_name)
        
        #if not threads:
        #    threads=self.threads
        #no threads in build index
        
        newOpts={"--":(fasta,),"-i":indexOut}
        mergedOpts={**newOpts,**kwargs}
        
        #call kallisto
        status=self.run_kallisto("index",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if status:
            #check if index file is present 
            if pu.check_files_exist(indexOut):
                self.kallisto_index=indexOut
                pu.print_green("kallisto_index is:"+self.kallisto_index)
                return True
        else:
            pu.print_boldred("Failed to create kallisto index")
            return False
    
    def perform_quant(self,sra_object,out_dir="",threads=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Run kallisto quant
        
        sra_object: SRA
            SRA object contatining paths to fastq files
        index_path: str
            path to the output directory
        index_name: str
            index name
        threads: int
            Number of threads
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to kallisto. This will override the existing options

        :return: Path to kallisto out directory
        :rtype: string
        """
        
        if not out_dir:
            out_dir=os.path.join(sra_object.location,"kallisto_out")
        
        if not threads:
            threads=self.threads
        
        
        if sra_object.layout == 'PAIRED':
            newOpts={"--threads":str(threads),"-o":out_dir,"--":(sra_object.localfastq1Path,sra_object.localfastq2Path),"-i":self.kallisto_index}
        else:
            newOpts={"--threads":str(threads),"-o":out_dir,"--single":"", "--":(sra_object.localfastqPath,),"-i":self.kallisto_index}
        
        
        #add input files to kwargs, overwrite newOpts if kwargs is present
        mergedOpts={**newOpts,**kwargs}
        
        #call kallisto
        status=self.run_kallisto("quant",verbose=verbose,quiet=quiet,logs=logs,objectid=sra_object.srr_accession,**mergedOpts)
        
        if status:
            #check if sam file is present in the location directory of sra_object
            if pu.check_files_exist(os.path.join(out_dir,"abundance.tsv")):
                return out_dir
        
        pu.print_boldred("kallisto quant failed")
        return ""
        
    
    
    def run_kallisto(self,subcommand,valid_args=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper for running kallisto.
        
        Parameters
        ----------
        
        subcommand: str
            subcommand for kallisto
        valid_args: list
            List of valid arguments, arguments in kwargs not in this list will be ignored
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to kallisto. This will override the existing options

        :return: Returns the status of kallisto. True is passed, False if failed.
        :rtype: bool
        """
        
        #check for a valid index
        if subcommand!="index":
            if not self.check_index():
                raise Exception("ERROR: Invalid kallisto index. Please run build index to generate an index.")
        
            
        kallisto_Cmd=['kallisto',subcommand]
        kallisto_Cmd.extend(pu.parse_unix_args(valid_args,kwargs))
        
        #start ececution
        status=pe.execute_command(kallisto_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,command_name=" ".join(kallisto_Cmd[0:2]))
        if not status:
            pu.print_boldred("kallisto failed")
        return status       
    
    def check_index(self):
        """Check valid kallisto index
        """
        if hasattr(self,'kallisto_index'):
            return(pu.check_files_exist(self.kallisto_index))
        return False
            


class Salmon(Quant):
    """This class represents salmon
    
    salmon_index: string
        Path to salmon index
    threads: int
        Number of threads
    """      
    def __init__(self,salmon_index,threads=None):    
        super().__init__() 
        self.programName="salmon"
        self.dep_list=[self.programName]        
        if not pe.check_dependencies(self.dep_list):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        """
        ##salmon index
        self.validArgsIndex=['-v','--version','-h','--help','-t','--transcripts','-k','--kmerLen','-i',
                             '--index','--gencode','--keepDuplicates','-p','--threads','--perfectHash',
                             '--type','-s','--sasamp']
        ##salmon quant read
        self.validArgsQuantReads=['--help-reads','-i','--index','-l','--libType','-r','--unmatedReads',
                                  '-1','--mates1','-2','--mates2','-o','--output','--discardOrphansQuasi',
                                  '--allowOrphansFMD','--seqBias','--gcBias','-p','--threads','--incompatPrior',
                                  '-g','--geneMap','-z','--writeMappings','--meta','--alternativeInitMode',
                                  '--auxDir','-c','--consistentHits','--dumpEq','-d','--dumpEqWeights',
                                  '--fasterMapping','--minAssignedFrags','--reduceGCMemory','--biasSpeedSamp',
                                  '--strictIntersect','--fldMax','--fldMean','--fldSD','-f','--forgettingFactor',
                                  '-m','--maxOcc','--initUniform','-w','--maxReadOcc','--noLengthCorrection',
                                  '--noEffectiveLengthCorrection','--noFragLengthDist','--noBiasLengthThreshold',
                                  '--numBiasSamples','--numAuxModelSamples','--numPreAuxModelSamples','--useVBOpt',
                                  '--rangeFactorizationBins','--numGibbsSamples','--numBootstraps','--thinningFactor',
                                  '-q','--perTranscriptPrior','--vbPrior','--writeOrphanLinks','--writeUnmappedNames',
                                  '-x','--quasiCoverage']
        ##salmon quant alignment
        self.validArgsQuantAlign=['--help-alignment','-l','--libType','-a','--alignments','-t','--targets','-p',
                                  '--threads','--seqBias','--gcBias','--incompatPrior','--useErrorModel',
                                  '-o','--output','--meta','-g','--geneMap','--alternativeInitMode','--auxDir'
                                  ,'--noBiasLengthThreshold','--dumpEq','-d','--dumpEqWeights','--fldMax',
                                  '--fldMean','--fldSD','-f','--forgettingFactor','--minAssignedFrags',
                                  '--gencode','--reduceGCMemory','--biasSpeedSamp','--mappingCacheMemoryLimit',
                                  '-w','--maxReadOcc','--noEffectiveLengthCorrection','--noFragLengthDist',
                                  '-v','--useVBOpt','--rangeFactorizationBins','--perTranscriptPrior','--vbPrior',
                                  '--numErrorBins','--numBiasSamples','--numPreAuxModelSamples','--numAuxModelSamples',
                                  '-s','--sampleOut','-u','--sampleUnaligned','--numGibbsSamples','--numBootstraps',
                                  '--thinningFactor']
        ##salmon quantmerge
        self.validArgsQuantMerge=['--quants','--names','-c','--column','-o','--output']

        self.valid_args=pu.get_union(self.validArgsIndex,self.validArgsQuantReads,self.validArgsQuantAlign,self.validArgsQuantMerge)
        """
        
        if not threads:
           threads=os.cpu_count()
         
        self.threads=threads
        
        #if index is passed, update the passed arguments
        if len(salmon_index)>0 and pu.check_salmon_index(salmon_index):
            print("salmon index is: "+salmon_index)
            self.salmon_index=salmon_index
        else:
            print("No salmon index provided. Please build index now to generate an index...")
            
            
            
    def build_index(self,index_path,index_name,fasta,threads=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """
        build salmon index and store the path to index in self
        
        index_path: str
            path to the output directory
        index_name: str
            index name
        fasta: str
            Path to fasta file
        threads: int
            Number of threads
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to salmon. This will override the existing options
            
        :return: status of salmon index
        :rtype: bool
        """
        
        #check input
        if not pu.check_files_exist(fasta):
            pu.print_boldred("{} does not exist. Exiting".format(fasta))
            return False
        #create out dir
        if not pu.check_paths_exist(index_path):
            if not pu.mkdir(index_path):
                print("ERROR in building hisat2 index. Failed to create index directory.")
                return False
        indexOut=os.path.join(index_path,index_name)
        
        if not threads:
            threads=self.threads
            
        newOpts={"--threads":str(threads),"-t":fasta,"-i":indexOut}
        
        mergedOpts={**newOpts,**kwargs}
        
        #call salmon
        status=self.run_salmon("index",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if status:
            #check if sam file is present in the location directory of sra_object
            #if check_files_exist(os.path.join(indexOut,"versionInfo.json")): #not sure if this is reliable
            if pu.check_paths_exist(indexOut):
                self.salmon_index=indexOut
                pu.print_green("salmon index is:"+self.salmon_index)
                return True
        
        pu.print_boldred("Failed to create salmon index")
        return False
        
        
    
    def perform_quant(self,sra_object,out_dir="",lib_type=None,threads=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """run salmon quant
        sra_object: SRA
            An SRA object with valid fastq files
        lib_type: str
            Library type. Default:A
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
            Options to pass to salmon. This will override the existing options

        :return: Path to salmon out directory
        :rtype: string
        """
        if not lib_type:
            lib_type="A"
            
        if not out_dir:
            out_dir=os.path.join(sra_object.location,"salmon_out")
        
        if not threads:
            threads=self.threads
        
        if sra_object.layout == 'PAIRED':
            newOpts={"--threads":str(threads),"-o":out_dir,"-l":lib_type,"-1":sra_object.localfastq1Path,"-2":sra_object.localfastq2Path,"-i":self.salmon_index}
        else:
            newOpts={"--threads":str(threads),"-o":out_dir,"-l":lib_type,"-r":sra_object.localfastqPath,"-i":self.salmon_index}
        
        
        #add input files to kwargs, overwrite newOpts with kwargs
        mergedOpts={**newOpts,**kwargs}
        
        #call salmon
        status=self.run_salmon("quant",verbose=verbose,quiet=quiet,logs=logs,objectid=sra_object.srr_accession,**mergedOpts)
        
        if status:
            #check if sam file is present in the location directory of sra_object
            if pu.check_files_exist(os.path.join(out_dir,"quant.sf")):
                return out_dir
        
        pu.print_boldred("salmon quant failed")
        return ""
        
        
        
    def run_salmon(self,subcommand,valid_args=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper for running salmon.
        
        Parameters
        ----------
        
        subcommand: str
            subcommand for salmon
        valid_args: list
            List of valid arguments
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to salmon. This will override the existing options

        :return: Returns the status of salmon. True is passed, False if failed.
        :rtype: bool
        """
        
        #check for a valid index
        if subcommand!="index":
            if not self.check_index():
                raise Exception("ERROR: Invalid salmon index. Please run build index to generate an index.")
        
                   
        salmon_Cmd=['salmon',subcommand]
        salmon_Cmd.extend(pu.parse_unix_args(valid_args,kwargs))
        
        #start ececution
        status=pe.execute_command(salmon_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,command_name=" ".join(salmon_Cmd[0:2]))
        if not status:
            pu.print_boldred("salmon failed")
        return status 

    def check_index(self):
        if hasattr(self,'salmon_index'):
            return pu.check_salmonindex(self.salmon_index)
        return False
    
    
