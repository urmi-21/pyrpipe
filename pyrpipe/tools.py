#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 14:54:22 2019

@author: usingh
"""
from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe
import os
import yaml

class RNASeqTools:
    def __init__(self):
        self.category="RNASeqTools"
    

class Samtools(RNASeqTools):
    """init function allows to specify commonly tuned parameters like threads, out_dir, memory etc.
    More specific parameters can be passed when calling the specific functions. There is option to override the parameters specified here in init later.
    
    threads: int
        Number of threads samtools will use
    max_memory: int
        Max memory to use in GB
    
    """
    def __init__(self,threads=None,max_memory=None):
        self.programName="samtools"
        #check if hisat2 exists
        if not pe.check_dependencies([self.programName]):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        self.threads=threads
        #Default: if threads are None use 80% of threads to avaoid memory issues
        if not self.threads:
            self.threads=int(os.cpu_count()*0.8)
            
        self.max_memory=max_memory
        
        
    def sam_to_bam(self,sam_file,out_dir="",out_suffix="",threads=None, delete_sam=False,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Convert sam file to a bam file. 
        Output bam file will have same name as input sam.
        sam_file: string
            Path to input Sam file
        out_suffix: string
            Suffix for the output sam file
        threads: int
            Number of threads. Default: Use self.threads initialized in init().
        delete_sam: bool
            delete the sam file after conversion
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to samtools. This will override the existing options 

        :return: Returns the path to the bam file. Returns empty string if operation failed.
        :rtype: string
        """        
        if not out_dir:            
            out_dir=pu.get_file_directory(sam_file)
        else:
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
        
        fname=pu.get_file_basename(sam_file)
        
        #output will be out_bam
        out_bam=os.path.join(out_dir,fname+out_suffix+'.bam')
        
        #handle threads
        if not threads:
            threads=self.threads
        
        newOpts={"--":(sam_file,),"-o":out_bam,"-@":str(threads),"-b":""}
        
        #add (and override) any arguments provided via kwargs
        mergedOpts={**newOpts,**kwargs}
        
        status=self.run_samtools("view",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
                
        if not status:
            print("Sam to bam failed for:"+sam_file)
            return ""
        
        #check if bam file exists
        if not pu.check_files_exist(out_bam):
            return ""
        
        #delete_sam_file
        if delete_sam:
            if not pe.deleteFileFromDisk(sam_file):
                print("Error deleting sam file:"+sam_file)
                
        #return path to file
        return out_bam
        
        
        
        
    #sort bam file.output will be bam_file_sorted.bam
    def sort_bam(self,bam_file,out_dir="",out_suffix="",threads=None,delete_bam=False,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Sorts an input bam file. Outpufile will end in _sorted.bam
        bam_file: str
            Path to the input bam file
        out_dir: str
            Path to output directory
        out_suffix: str
            Output file suffix
        threads: int
            Number of threads. Default: Use self.threads initialized in init().
        delete_bam: bool
            Delete input bam_file
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to samtools. This will override the existing options 

        :return: Returns path to the sorted bam file. Returns empty string if operation failed.
        :rtype: string
        
        """
        if not out_dir:
            out_dir=pu.get_file_directory(bam_file)
        else:
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
                
        fname=pu.get_file_basename(bam_file)
        #output will be out_bam
        outSortedbam_file=os.path.join(out_dir,fname+out_suffix+'_sorted.bam')
        
        #handle threads
        if not threads:
            threads=self.threads
        
        newOpts={"--":(bam_file,),"-o":outSortedbam_file,"-@":str(threads)}
        mergedOpts={**newOpts,**kwargs}
        
        status=self.run_samtools("sort",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if not status:
            print("Bam sort failed for:"+bam_file)
            return ""
        
        #check if bam file exists
        if not pu.check_files_exist(outSortedbam_file):
            return ""

        if delete_bam:
            if not pe.deleteFileFromDisk(bam_file):
                print("Error deleting sam file:"+bam_file)
                
        #return path to file
        return outSortedbam_file
    
    def sam_sorted_bam(self,sam_file,out_dir="",out_suffix="",threads=None,delete_sam=False,delete_bam=False,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Convert sam file to bam and sort the bam file.
        sam_file: str
            Path to the input sam file
        out_dir: str
            Path to output directory
        out_suffix: str
            Output file suffix
        threads: int
            Number of threads. Default: Use self.threads initialized in init().
        delete_sam: bool
            Delete input sam_file
        delete_bam: bool
            Delete the intermediate unsorted bam_file
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to samtools. This will override the existing options 

        :return: Returns path to the sorted bam file. Returns empty string if operation failed.
        :rtype: string
        """
        
        sam2bam_file=self.sam_to_bam(sam_file,threads=threads,delete_sam=delete_sam,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**kwargs)
        
        if not sam2bam_file:
            return ""
            

        bamSorted=self.sort_bam(sam2bam_file,out_dir, out_suffix,threads=threads,delete_bam=delete_bam,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**kwargs)
        
        if not bamSorted:
            return ""
        
        return bamSorted
    
    
    def merge_bam(self,*args,out_file="merged",out_dir="",threads=None,delete_bams=False,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Merge multiple bam files into a single file
        
        Parameters
        ----------
        args: *args
            Input bam files to merge
        out_file: string
            Output file name to save the results. .bam will be added at the end.
        out_dir: string
            Path where to save the merged bam file. Default path is the same as the first bam_file's
        threads: int
            Number of threads. Default: Use self.threads initialized in init().
        delete_bams: bool
            Delete input bam files after merging.
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to samtools. This will override the existing options 

        :return: Returns the path to the merged bam file.
        :rtype: string
        """
               
        if len(args) < 2:
            print("Please supply at least 2 files to merge")
            return ""
        
        if not out_dir:
            out_dir=pu.get_file_directory(args[0])
        else:
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
        
        outMergedFile=os.path.join(out_dir,out_file+".bam")
        
        #handle threads
        if not threads:
            threads=self.threads
        
        newOpts={"-@":str(threads),"--":(outMergedFile,)+args}
        
        #override parameters by supplying as kwargs
        mergedOpts={**newOpts,**kwargs}
        
        status=self.run_samtools("merge",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if not status:
            print("Bam merge failed for:"+outMergedFile)
            return ""
        
        #check if bam file exists
        if not pu.check_files_exist(outMergedFile):
            return ""
        

        if delete_bams:
            for bam_file in args:
                if not pe.deleteFileFromDisk(bam_file):
                    print("Error deleting sam file:"+bam_file)
                    
        return outMergedFile
        
        
        
    def run_samtools(self,sub_command,valid_args=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """A wrapper to run samtools.
        
        Parameters
        ----------
        
        sub_command: string
            sub_command to pass to samtools e.g. sort, merge etc
        valid_args: list
            A list containing valid parameters. Parameters in kwargs not in this list will be ignored. Default: None
        arg1: dict
            arguments to pass to samtools. 
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to samtools. This will override the existing options 

        :return: Returns the status of samtools. True is passed, False if failed.
        :rtype: bool
        """
            
       
       
        samtools_cmd=['samtools',sub_command]
        #add options
        samtools_cmd.extend(pu.parse_unix_args(valid_args,kwargs))
                
        #start ececution
        status=pe.execute_command(samtools_cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            pu.print_boldred("samtools failed")
        
        #return status
        return status
        
        
        
        
        
class Portcullis(RNASeqTools):
    def __init__(self,threads=None,max_memory=None):
        self.programName="portcullis"
        self.dep_list=[self.programName]
        #check if program exists
        if not pe.check_dependencies(self.dep_list):
            raise Exception("ERROR: "+ self.programName+" not found.")
               
        #use max threads by default
        if not threads:
            threads=os.cpu_count()
        self.threads=threads
            
        
        
    def run_portcullisFull(self,reference_fasta,bam_file,threads=None,out_dir="",delete_bam=False,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """
        run portculis full
        
        Parameters
        ----------
        
        reference_fasta: string
            Path to the reference fasta file
        bam_file: string
            Path to input bam file
        threads: int
            Number of threads
        out_dir: string
            Path to the out put dir. current directory is not given.
        delete_bam: bool
            delete input bam file        
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to portcullis. This will override the existing options 
        
        """
        
        if not pu.check_files_exist(reference_fasta,bam_file):
            print ("Please check input for portcullis.")
            return ""
        
        #handle threads
        if not threads:
            threads=self.threads
            
        #add out dir path
        if not out_dir:
            out_dir=os.path.join(os.getcwd(),"portcullis_out")
            
        newOpts={"--":(reference_fasta,bam_file),"--threads":str(threads),"-o":out_dir}
        
        mergedOpts={**newOpts,**kwargs}
                
        status=self.run_portcullis("full",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if not status:
            print("portcullis full failed for:"+bam_file)
            return ""
        
        #check if bam file exists
        if not pu.check_paths_exist(out_dir):
            return ""

        if delete_bam:
            if not pe.deleteFileFromDisk(bam_file):
                    print("Error deleting bam file:"+bam_file)
        
        return out_dir
    
    def run_portcullis(self,sub_command,valid_args=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """
        Wrapper to run portcullis.
        
        Parameters
        ----------
        
        sub_command: string
            sub_command to pass to portcullis e.g. full, prep, junc etc.
        valid_args: list
            A list of valid arguments. Arguments outside this list will be ignored. If empty or None, accepts all arguments.
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            arguments to pass to portcullis. 

        :return: Returns the status of portcullis. True is passed, False if failed.
        :rtype: bool
        """
        
        
        portcullis_cmd=['portcullis',sub_command]
        #add options
        portcullis_cmd.extend(pu.parse_unix_args(valid_args,kwargs))

        #start ececution
        status=pe.execute_command(portcullis_cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            pu.print_boldred("portcullis failed")
                
        #return status
        return status
        
        
        

class Mikado(RNASeqTools):
    """Mikado constructor
    threads: int
        Number of threads to use
    max_memory: float
        Max memory to use in (GB)
    """
    def __init__(self,threads=None,max_memory=None):
        self.programName="mikado"
        self.dep_list=[self.programName]
        #check if program exists
        if not pe.check_dependencies(self.dep_list):
            raise Exception("ERROR: "+ self.programName+" not found.")
    
       
        #use max threads by default
        if not threads:
            threads=os.cpu_count()
            
        self.threads=threads
        
        
 
    def createMikadoGTFlist(self,out_file,out_dir,searchPath,searchQuery="*.gtf",strand=False):
        """Create a file to be used by mikado configure
        out_file: str
            outfile name
        out_dir: str
            path to out_dir
        searchPath: str
            Path where gtf/gff files will be searched
        searchQuery: str
            Query to perform search. Default: "*.gtf"
        strand: bool
            Stranded flag: Default false
        
            
        """
        
        files=pe.find_files(searchPath,searchQuery,recursive=True)
        args=files
        
        #create out dir
        if not pu.check_paths_exist(out_dir):
            pu.mkdir(out_dir)
        outFilePath=os.path.join(out_dir,out_file+".txt")
        
        
        gtfs=[]
        for l in args:
            thisName=pu.get_file_basename(l)
            if thisName:
                gtfs.append("\t".join([l,thisName,str(strand)]))
        
        f=open(outFilePath,"w")
        f.write("\n".join(gtfs))
        f.close()
        
        pu.print_green("Mikado list file written to:"+outFilePath)
        return outFilePath
                

        
    def runMikadoFull(self,listFile,genome,mode,scoring,junctions,config_out_file,blast_targets,blastx_object,out_dir=None,threads=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Run whole mikado pipeline
        Output will be stored to out_dir/
        
        listFile: str
            Input list file to mikado configure
        genome: str
            Path to gename fasta file
        mode: str
            mikado mode
        scoring: str
            scoring file
        junctions: str
            path to junctions file
        config_out_file:str
            path to output configure file
        balst_targets: str
            Path to blast targets
        blastx_object: object
            A pyrpipe object to use homology search e.g. Diamond
        out_dir: str
            path to out directory
        threads: int
            Number of threads to use
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            arguments to pass to mikado. 
            
        """
        
        #use threads if provided
        if not threads:
            threads=self.threads
        
        #run mikado config
        config_file=self.runMikadoConfigure(listFile,genome,mode,scoring,junctions,config_out_file, threads=threads, out_dir=out_dir,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**kwargs)
        if not pu.check_files_exist(config_file):
            pu.print_boldred("Mikado configure failed")
            return False
        
        #run mikado prep
        self.runMikadoPrepare(config_file,threads=threads, out_dir=out_dir,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**kwargs)
        
        
        mikado_prep_fa=os.path.join(out_dir,"mikado_prepared.fasta")
        
        #run blast to get blast/diamond xml
        
        xml=blastx_object.run_align(mikado_prep_fa,"mikado_prep_blast.xml",command="blastx",out_fmt=5,fmt_string=None,out_dir=out_dir,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        
        #print("XML:"+xml)
        
        #run transdecoder/prodigal to get orfs
        
        txd=Transdecoder()
        longOrfOut=txd.run_transdecoder_longorfs(mikado_prep_fa,out_dir=out_dir+"/longorfsout")
        preddir=out_dir+"/predout"
        #predout=txd.run_transdecoder_predict(mikado_prep_fa,longOrfOut,out_dir=preddir)
        #print(predout)
        
        orfs=os.path.join(preddir,"mikado_prepared.fasta.transdecoder.bed")
        
        #run mikado ser
        self.runMikadoSerialise(config_file,blast_targets,orfs, xml,threads=threads,out_dir=out_dir,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**kwargs)
        
        #run mikado pick
        self.runMikadoPick(config_file,threads=threads,out_dir=out_dir,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**kwargs)
        
        #return dir containing mikado results
        return out_dir
        
    
    
    
    def runMikadoConfigure(self,listFile,genome,mode,scoring,junctions,out_file,threads=None,out_dir=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper to run mikado configure
        Make sure the paths in list file are global.
        
        Parameters
        ----------

        :return: Path to the created configuration file
        :rtype: string
        """
        
        #check all file exists
        if not pu.check_files_exist(listFile,genome,junctions):
            print("Please check mikado input")
            return ""
        
        #create out dir
        if out_dir==None:
            out_dir=os.getcwd()
        if not pu.check_paths_exist(out_dir):
            if not pu.mkdir(out_dir):
                raise Exception("Exception in mikado configure.")
            
        outFilePath=os.path.join(out_dir,out_file+".yaml")
        
        if not threads:
            threads=self.threads
        
        newOpts={"--threads":str(threads),"--list":listFile,"--reference":genome,"--mode":mode,"--scoring":scoring,"--junctions":junctions,"--out-dir":out_dir,"--":(outFilePath,)}
        
        #merge with kwargs
        mergedOpts={**newOpts,**kwargs}
        
        status=self.runMikado("configure",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if not status:
            pu.print_boldred("Mikado configure failed.\nPlease make sure the paths in list file are global.")
            return ""
        
        #check if bam file exists
        if not pu.check_files_exist(outFilePath):
            return ""
        
        
        return outFilePath
        
    
    def runMikadoPrepare(self,yamlconf, threads=None,out_dir=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper to run mikado prepare
        """
        
        #check input files exist
        if not pu.check_files_exist(yamlconf):
            print("Please check the input configuration to mikado.")
            return ""
        if not out_dir:
            out_dir=os.getcwd()
            
        if not threads:
            threads=self.threads

        newOpts={"--procs":str(threads),"--output-dir":out_dir,"--json-conf":yamlconf}
        
        #merge with kwargs
        mergedOpts={**newOpts,**kwargs}
        
        status=self.runMikado("prepare",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if not status:
            print("Mikado prepare failed for:"+yamlconf)
            return ""
        
        #check if bam file exists
        if not pu.check_paths_exist(out_dir):
            return ""
        
        #after running mikado prep update the mikado config file to contain absolute path of mikafo_prepared.fa
        #otherwise serialize step will fail if out dirs are different
        #read the cofig.json file
        with open(yamlconf) as f:
            conf_data = yaml.safe_load(f)
        tx=conf_data['serialise']['files']['transcripts']
        conf_data['serialise']['files']['transcripts']=os.path.join(out_dir,tx)
        #update input to pick
        gtf=conf_data['pick']['files']['input']
        conf_data['pick']['files']['input']=os.path.join(out_dir,gtf)
        #write to file
        with open(yamlconf,'w') as f:
            yaml.dump(conf_data, f)
       
        
        return out_dir
        
        
        
    def runMikadoSerialise(self,yamlconf,blast_targets,orfs,xml,threads=None,out_dir=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper to run mikado serialise
        """
        #check input files exist
        if not pu.check_files_exist(blast_targets,orfs,xml):
            print("Please check the input to mikado.")
            return ""
        if not out_dir:
            out_dir=os.getcwd()
        
        if not threads:
            threads=self.threads
        
        newOpts={"--procs":str(threads),"--json-conf":yamlconf,"--blast_targets":blast_targets,"--xml":xml,"--orfs":orfs,"--output-dir":out_dir}
        
        #merge with kwargs
        mergedOpts={**newOpts,**kwargs}
        
        
        
        status=self.runMikado("serialise",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if not status:
            print("Mikado serialise failed for:"+yamlconf)
            return ""
        
        #check if out path exists
        if not pu.check_paths_exist(out_dir):
            return ""
        
        #after running mikado serialize update the mikado config file to contain absolute path of mikadodb
        #otherwise pick step will fail
        #read the cofig.json file
        with open(yamlconf) as f:
            conf_data = yaml.safe_load(f)
        db=conf_data['db_settings']['db']
        conf_data['db_settings']['db']=os.path.join(out_dir,pu.get_filename(db))
        
        
        #write to file
        with open(yamlconf,'w') as f:
            yaml.dump(conf_data, f)
        
        return out_dir
        
        
    def runMikadoPick(self,yamlconf,threads=None,out_dir=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper to run mikado pick
        """
        #check input files exist
        if not pu.check_files_exist(yamlconf):
            print("Please check the input to mikado.")
            return ""
        if not out_dir:
            out_dir=os.getcwd()
            
        if not threads:
            threads=self.threads
        
        newOpts={"--procs":str(threads),"--json-conf":yamlconf,"--output-dir":out_dir}
        
        #merge with kwargs
        mergedOpts={**newOpts,**kwargs}
        
        status=self.runMikado("pick",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if not status:
            print("Mikado pick failed for:"+yamlconf)
            return ""
        
        #check if bam file exists
        if not pu.check_paths_exist(out_dir):
            return ""
        
        return out_dir
        
        
    def runMikado(self,sub_command,valid_args=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper to run mikado
        """
        valid_commands=['configure','prepare','serialise','pick','compare']
        if sub_command not in valid_commands:
            pu.print_boldred("Invalid command: "+sub_command+". Exiting...")
            return False
        
      
        mikado_Cmd=['mikado',sub_command]
        #add options
        mikado_Cmd.extend(pu.parse_unix_args(valid_args,kwargs))
                
        #print("Executing:"+" ".join(mergedArgsDict))
        
        #start ececution
        status=pe.execute_command(mikado_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            pu.print_boldred("mikado failed")
        #return status
        return status
        
        
        
        
        
        
class Ribocode(RNASeqTools):
    def __init__(self,threads=None):
        self.programName="RiboCode"
        self.dep_list=[self.programName]
        #check if program exists
        if not pe.check_dependencies(self.dep_list):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        #use max threads by default
        if not threads:
            threads=os.cpu_count()
        self.threads=threads
        
        
        
    def runRibocode(self,gtf,genome,bam,l="no",outsuffix="ribocode_out",verbose=False,quiet=False,logs=True,objectid="NA"):
        """Wrapper to run ribocode in one step
        """
        
        #check input
        if not pu.check_files_exist(gtf,genome,bam):
            pu.print_boldred("Please check input files for Ribocode")
            return ""
        
        out_dir=pu.get_file_directory(gtf)
        outFile=os.path.join(out_dir,outsuffix)
        
        newOpts={"-g":gtf,"f":genome,"-r":bam,"-l":l,"-o":outFile}
        
        ribocode_Cmd=['RiboCode_onestep']
        ribocode_Cmd.extend(pu.parse_unix_args(self.valid_args,newOpts))
        
        status=pe.execute_command(ribocode_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            pu.print_boldred("ribocode failed")
            return ""
        
        return outFile
        
        
        

        
       
        
        
        
class Diamond(RNASeqTools):
    def __init__(self,index,threads=None,mode=None):
        self.programName="diamond"
        self.dep_list=[self.programName]
        #check if program exists
        if not pe.check_dependencies(self.dep_list):
            raise Exception("ERROR: "+ self.programName+" not found.")
        """
        self.valid_args=['-p','--db','-d','--out','-o','--outfmt','-f','--verbose','--log','--quiet','--in','--query','-q','--strand','--un','--al','--unal','--max-target-seqs','-k','--top','--range-culling','--compress','--evalue','-e',
                         '--min-score','--id','--query-cover','--subject-cover','--sensitive','--more-sensitive','--block-size','-b','--index-chunks','-c','--tmpdir','-t','--gapopen','--gapextend','--frameshift','-F','--long-reads','--matrix','--custom-matrix',
                         '--lambda','--K','--comp-based-stats','--masking','--query-gencode','--salltitles','--sallseqid','--no-self-hits','--taxonmap','--taxonnodes','--taxonlist','--algo','--bin','--min-orf','-l','--freq-sd','--id2','--window','-w',
                         '--xdrop','-x','--ungapped-score','--hit-band','--hit-score','--gapped-xdrop','-X','--band','--shapes','-s','--shape-mask','--index-mode','--rank-ratio','--rank-ratio2','--max-hsps','--range-cover','--dbsize','--no-auto-append',
                         '--xml-blord-format','--daa','-a','--forwardonly','--seq']     
        """
        self.valid_commands=['makedb','blastp','blastx','view','help','version','getseq','dbinfo']
        
        
        
        #use max threads by default
        if not threads:
            threads=os.cpu_count()    
        self.threads=threads
            
        #select mode
        valid_modes=['fast','sensitive','more-sensitive']
                    
        if mode in valid_modes:
            mode='--'+mode
        else:
            mode='--fast'
            
        self.mode=mode
        
        #check index
        self.index=index
        if not self.check_index():
            print("No valid index provided. Please build index...")
        else:
            self.index=index
            
    def build_index(self,in_fasta,dbname,out_dir=None,threads=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Build a diamond index and store its path in self
        """
        
        #check input files
        if not pu.check_files_exist(in_fasta):
            pu.print_boldred("Input fasta: {} not found...\n diamond makedb failed".format(in_fasta))
            return False
        #create out_dir
        if not out_dir:
            out_dir=os.getcwd()
        if not pu.check_paths_exist(out_dir):
            pu.mkdir(out_dir)
        
            
        #check if index already exists
        index_path=os.path.join(out_dir,dbname)
        self.index=index_path
        if self.check_index():
            pu.print_green("Diamond index: {} exists, using it...".format(self.index))
            self.index=index_path
            return True
        
        if not threads:
            threads=self.threads
        
        newOpts={"--in": in_fasta, "-d": index_path, "--threads":str(threads)}
        
        #add input files to kwargs, overwrite newOpts with kwargs 
        mergedOpts={**newOpts,**kwargs}
        
        #call run_diamond
        status=self.run_diamond("makedb",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if status:
            self.index=index_path
            return True
        
        return False
    
    
    def run_align(self,query,out_file,command=None,out_fmt=None,fmt_string=None,out_dir=None,threads=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """
        Parameters
        ----------
        query: string
            Path to query fasta file
        out_file: string
            name of output file
        command: string
            "blastx" or "blastp"
        out_fmt: int
            output format as specified by diamond: 
                0 = BLAST pairwise
                5 = BLAST XML
                6 = BLAST tabular
                100 = DIAMOND alignment archive (DAA)
                101 = SAM
        fmt_string: string
             space-separated list of keywords for out format 6
        out_dir: string
            path to out_directory
            Default: current working dir
        threads: int
            number of threads to use 
            Default:2
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        
        kwargs: dict
            arguments to pass to hisat2. This will override parametrs already existing in the self.passedArgumentList list but NOT replace them.

        :return: Returns the status of diamond. True is passed, False if failed.
        :rtype: bool
        """
        if not command:
            command="blastx"
        if command not in ['blastx','blastp']:
            pu.print_boldred("Invalid command: {}. Exiting...".format(command))
            return ""
        
        #check input
        if not pu.check_files_exist(query):
            pu.print_boldred("Input fasta query file: {} not found...\n diamond blastx failed".format(query))
            return ""
        if not out_dir:
            out_dir=os.getcwd()
        if not pu.check_paths_exist(out_dir):
            pu.mkdir(out_dir)
            
        if not out_fmt:
            out_fmt=0
        if out_fmt ==6:
            if not fmt_string:
                fmt_string=""
            out_fmt=str(out_fmt)+" "+fmt_string
        
        
        out_file_path=os.path.join(out_dir,out_file)
        
        if not threads:
            threads=self.threads
        
        newOpts={"-d": self.index, "-q": query, "-o":out_file_path,"-f":str(out_fmt),"--threads":str(threads)}
        
        #add input files to kwargs, overwrite kwargs with newOpts
        mergedOpts={**newOpts,**kwargs}
        
        #call run
        status=self.run_diamond(command,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if status:
            #check if sam file is present in the location directory of sra_object
            if pu.check_files_exist(out_file_path):
                return out_file_path
        else:
            return ""
        
        
        
    def run_diamond(self,subcommand,valid_args=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper for running diamond.
        
        Parameters
        ----------
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        
        kwargs: dict
            arguments to pass to hisat2. This will override parametrs already existing in the self.passedArgumentList list but NOT replace them.

        :return: Returns the status of diamond. True is passed, False if failed.
        :rtype: bool
        """
        
        #check for a valid index
        if subcommand=="blastx" or subcommand=="blastp":
            if not self.check_index():
                raise Exception("ERROR: Invalid Diamond index. Please run build_index() to generate an index.")
            
              
        diamond_cmd=['diamond',subcommand]
        #add options
        diamond_cmd.extend(pu.parse_unix_args(valid_args,kwargs))        
        
        #execute command
        cmd_status=pe.execute_command(diamond_cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not cmd_status:
            print("Diamond failed:"+" ".join(diamond_cmd))
     
        #return status
        return cmd_status
    
        
    def check_index(self):
        """Check a diamond index
        """
        if not hasattr(self,"index"):
            return False
        
        if pu.check_files_exist(self.index):
            return True
        
        if pu.check_files_exist(self.index+".dmnd"):
            return True
        
        return False
        
        
class Transdecoder(RNASeqTools):
    def __init__(self):
        self.programName="TransDecoder.LongOrfs"
        self.dep_list=["TransDecoder.LongOrfs" ,"TransDecoder.Predict"]
        #check if program exists
        if not pe.check_dependencies(self.dep_list):
            raise Exception("ERROR: "+ str(self.dep_list)+" not found.")
        
        
    def run_transdecoder_longorfs(self,infasta,out_dir=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        
        if not pu.check_files_exist(infasta):
            pu.print_boldred("Please check input file:"+infasta)
        
        
        if not out_dir:
            out_dir=os.getcwd()
            
            
        newOpts={"-t":infasta,"-O":out_dir}
        mergedOpts={**newOpts,**kwargs}
        
        #execute LongOrfs
        status=self.run_transdecoder('TransDecoder.LongOrfs',verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        if not status:
            pu.print_boldred("Transdecoder failed")
            return ""
   
        return out_dir
        
        
    def run_transdecoder_predict(self,infasta,longorfs_dir,out_dir=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        
        if not pu.check_files_exist(infasta):
            pu.print_boldred("Please check input file:"+infasta)
        if not pu.check_paths_exist(longorfs_dir):
            pu.print_boldred("Path {} doesn't exist".format(longorfs_dir))
            
        move_flag=True
        if not out_dir:
            out_dir=os.getcwd()
            move_flag=False
            
        if not pu.check_paths_exist(out_dir):
            pu.mkdir(out_dir)
            
        newOpts={"-t":infasta,"-O":longorfs_dir}
        mergedOpts={**newOpts,**kwargs}
        
        #execute Predict
        status=self.run_transdecoder('TransDecoder.Predict',verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        if not status:
            pu.print_boldred("Transdecoder failed")
            return ""
        
        #move output files to outdir
        if move_flag:
            outfile_prefix=pu.get_filename(infasta)+".transdecoder"
            pe.move_file(outfile_prefix+".bed",os.path.join(out_dir,outfile_prefix+".bed"),verbose)
            pe.move_file(outfile_prefix+".cds",os.path.join(out_dir,outfile_prefix+".cds"),verbose)
            pe.move_file(outfile_prefix+".gff3",os.path.join(out_dir,outfile_prefix+".gff3"),verbose)
            pe.move_file(outfile_prefix+".pep",os.path.join(out_dir,outfile_prefix+".pep"),verbose)
        return out_dir
        
        
    
    def run_transdecoder(self,command,valid_args=None,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper for running transdecoder.
        
        Parameters
        ----------
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        
        kwargs: dict
            arguments to pass to hisat2. This will override parametrs already existing in the self.passedArgumentList list but NOT replace them.

        :return: Returns the status of diamond. True is passed, False if failed.
        :rtype: bool
        """
        
        txd_cmd=[command]
        #add options
        txd_cmd.extend(pu.parse_unix_args(valid_args,kwargs))        
        
        #execute command
        cmd_status=pe.execute_command(txd_cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not cmd_status:
            print("Transdecoder failed:"+" ".join(txd_cmd))
     
        #return status
        return cmd_status


