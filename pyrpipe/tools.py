#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 14:54:22 2019

@author: usingh
"""
from pyrpipe.pyrpipe_utils import *
from pyrpipe.pyrpipe_engine import *

class RNASeqTools:
    def __init__(self):
        self.category="RNASeqTools"
    

##Replace with pysam

class Samtools(RNASeqTools):
    def __init__(self,**kwargs):
        self.programName="samtools"
        #check if hisat2 exists
        if not check_dependencies([self.programName]):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        self.validArgsList=['-b','-C','-1','-u','-h','-H','-c','-o','-U','-t','-L','-r',
                            '-R','-q','-l','-m','-f','-F','-G','-s','-M','-x','-B','-?','-S','-O','-T','-@']
        
        self.passedArgumentDict=kwargs
        
        
        
    def samToBam(self,samFile,outFileSuffix="",deleteSam=False,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
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
        
        status=self.runSamtools("view",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
                
        if not status:
            print("Sam to bam failed for:"+samFile)
            return ""
        
        #check if bam file exists
        if not check_files_exist(outBamFile):
            return ""
        
        #deletesamfile
        if deleteSam:
            if not deleteFileFromDisk(samFile):
                print("Error deleting sam file:"+samFile)
                
        #return path to file
        return outBamFile
        
        
        
        
    #sort bam file.output will be bamFile_sorted.bam
    def sortBam(self,bamFile,outFileSuffix="",deleteOriginalBam=False,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
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
        
        status=self.runSamtools("sort",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if not status:
            print("Bam sort failed for:"+bamFile)
            return ""
        
        #check if bam file exists
        if not check_files_exist(outSortedBamFile):
            return ""

        if deleteOriginalBam:
            if not deleteFileFromDisk(bamFile):
                print("Error deleting sam file:"+bamFile)
                
        #return path to file
        return outSortedBamFile
    
    def samToSortedBam(self,samFile,outFileSuffix="",deleteSam=False,deleteOriginalBam=False,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Convert sam file to bam and sort the bam file.
        
        Returns
        -------
        string
                Returns path to the sorted bam file. Returns empty string if operation failed.
        """
        
        sam2BamFile=self.samToBam(samFile,deleteSam=deleteSam,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**kwargs)
        
        if not sam2BamFile:
            return ""
            

        bamSorted=self.sortBam(sam2BamFile,outFileSuffix,deleteOriginalBam,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**kwargs)
        
        if not bamSorted:
            return ""
        
        return bamSorted
    
    
    def mergeBamFiles(self,*args,outFileName="merged",outPath="",deleteOriginalBamFiles=False,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Merge multiple bam files into a single file
        
        Parameters
        ----------
        outFileName: string
            Output file name to save the results. .bam will be added at the end.
        args:tuple
            Paths to bam files to combine
        outPath: string
            Path where to save the merged bam file. Default path is the same as the first bamfile's
        kwargs: dict
            arguments passed to samtools merge command
            
        Returns
        -------
        string
            Returns the path to the merged bam file.
        """
        
        if len(args) < 2:
            print("Please supply at least 2 files to merge")
            return ""
        
        if outPath == "":
            outPath=getFileDirectory(args[0])
        
        outMergedFile=os.path.join(outPath,outFileName+".bam")
        
        newOpts={"--":(outMergedFile,)+args}
        
        mergedOpts={**kwargs,**newOpts}
        
        status=self.runSamtools("merge",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if not status:
            print("Bam merge failed for:"+outMergedFile)
            return ""
        
        #check if bam file exists
        if not check_files_exist(outMergedFile):
            return ""
        

        if deleteOriginalBamFiles:
            for bamFile in args:
                if not deleteFileFromDisk(bamFile):
                    print("Error deleting sam file:"+bamFile)
                    
        return outMergedFile
        
        
        
    def runSamtools(self,subCommand,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """A wrapper to run samtools.
        
        Parameters
        ----------
        subcommand: string
            Subcommand to pass to samtools e.g. sort, merge etc
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
        samtools_Cmd.extend(parse_unix_args(self.validArgsList,mergedArgsDict))
                
        #start ececution
        status=executeCommand(samtools_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            print_boldred("samtools failed")
        
        #return status
        return status
        
        
        
        
        
class Portcullis(RNASeqTools):
    def __init__(self,**kwargs):
        self.programName="portcullis"
        self.depList=[self.programName]
        #check if program exists
        if not check_dependencies(self.depList):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        self.validArgsList=['-t','--threads','-v','--verbose','--help','-o','-b',
                            '--bam_filter','--exon_gff','--intron_gff','--source',
                            '--force','--copy','--use_csi','--orientation','--strandedness',
                            '--separate','--extra','-r','--max_length','--canonical','--min_cov',
                            '--save_bad']
        
        self.passedArgumentDict=kwargs
        
        
    def runPortcullisFull(self,referenceFasta,bamFile,outDir="",deleteOriginalBamFile=False,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """
        run portculis full
        
        Parameters
        ----------
        referenceFasta: string
            Path to the reference fasta file
        bamFile: string
            Path to input bam file
        outDir: string
            Path to the out put dir. current directory is not given.
        """
        
        if not check_files_exist(referenceFasta,bamFile):
            print ("Please check input for portcullis.")
            return ""
        
        
        newOpts={"--":(referenceFasta,bamFile)}
        mergedOpts={**kwargs,**newOpts}
        #add out dir path
        if not outDir:
            outDir=os.path.join(os.getcwd(),"portcullis_out")
                  
        mergedOpts={**mergedOpts,**{"-o":outDir}}
        
        status=self.runPortcullis("full",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if not status:
            print("portcullis full failed for:"+bamFile)
            return ""
        
        #check if bam file exists
        if not checkPathsExists(outDir):
            return ""

        if deleteOriginalBamFile:
            if not deleteFileFromDisk(bamFile):
                    print("Error deleting bam file:"+bamFile)
        
        return outDir
    
    def runPortcullis(self,subCommand,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """
        Wrapper to run portcullis.
        
        Parameters
        ----------
         subcommand: string
            Subcommand to pass to portcullis e.g. full, prep, junc etc.
        arg1: dict
            arguments to pass to samtools. This will override parametrs already existing in the self.passedArgumentDict list but NOT replace them.
            
        Returns
        -------
        bool:
                Returns the status of portcullis. True is passed, False if failed.
        """
        
        
        #override existing arguments
        mergedArgsDict={**self.passedArgumentDict,**kwargs}
       
        portcullis_Cmd=['portcullis',subCommand]
        #add options
        portcullis_Cmd.extend(parse_unix_args(self.validArgsList,mergedArgsDict))
                
        print("Executing:"+" ".join(portcullis_Cmd))
        
        
        #start ececution
        status=executeCommand(portcullis_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            print_boldred("portcullis failed")
                
        #return status
        return status
        
        
        

class Mikado(RNASeqTools):
    def __init__(self,**kwargs):
        self.programName="mikado"
        self.depList=[self.programName]
        #check if program exists
        if not check_dependencies(self.depList):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        self.validArgsList=[]        
        self.passedArgumentDict=kwargs
        
        
    #mikado prepare --start-method spawn --log mikadoPrepareLog --minimum_length 30
    #--procs 28 --output-dir mikadoPrepOutNew 
    #--fasta /pylon5/mc5pl7p/usingh/lib/hisatIndex/ensembl_release98/Homo_sapiens.GRCh38.dna.primary_assembly.fa 
    #--list smolGtfList
    
    def searchGTFtolist(self, outFileName,searchPath=os.getcwd(),searchQuery="*.gtf",outDir=os.getcwd(),strand=False):
        searchCmd=['find',searchPath,'-name',searchQuery]
        st=runLinuxCommand(searchCmd)
        if st[0]==0:
            output=st[1].decode("utf-8").split("\n")
            
        return self.createMikadoGTFlist(outFileName,*output,outDir=outDir,strand=strand)


        
    
    def createMikadoGTFlist(self,outFileName,*args,outDir=os.getcwd(),strand=False):
        """Create a file to be used by mikado configure
        """
        
        outFilePath=os.path.join(outDir,outFileName+".yaml")
        
        
        gtfs=[]
        for l in args:
            thisName=getFileBaseName(l)
            if thisName:
			    #print("\t".join([l,thisName,strand]))
                gtfs.append("\t".join([l,thisName,str(strand)]))
        
        f=open(outFilePath,"w")
        f.write("\n".join(gtfs))
        f.close()
        
        print_green("Mikado list file written to:"+outFilePath)
        return outFilePath
                

        
    def runMikadoFull(self):
        """Run whole mikado pipeline
        """
        pass
    
    def runMikadoConfigure(self,listFile,genome,mode,scoring,junctions,outFileName,outDir=os.getcwd(),verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper to run mikado configure
        Make sure the paths in list file are global.
        Parameters
        ----------
        
        returns
        -------
        string
            Path to the created configuration file
        """
        
        #check all file exists
        if not check_files_exist(listFile,genome,junctions,scoring):
            print("Please check mikado input")
            return ""
        
        #create out dir
        if not checkPathsExists(outDir):
            if not mkdir(outDir):
                raise Exception("Exception in mikado configure.")
            
        outFilePath=os.path.join(outDir,outFileName+".yaml")
        
        newOpts={"--list":listFile,"--reference":genome,"--mode":mode,"--scoring":scoring,"--junctions":junctions,"--":(outFilePath,)}
        
        #merge with kwargs
        mergedOpts={**kwargs,**newOpts}
        
        status=self.runMikado("configure",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if not status:
            print_boldred("Mikado configure failed.\nPlease make sure the paths in list file are global.")
            return ""
        
        #check if bam file exists
        if not check_files_exist(outFilePath):
            return ""
        
        return outFilePath
        
    
    def runMikadoPrepare(self,jsonconf, outDir="",verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper to run mikado prepare
        """
        
        #check input files exist
        if not check_files_exist(jsonconf):
            print("Please check the input configuration to mikado.")
            return ""
        if not outDir:
            outDir=os.getcwd()

        newOpts={"--output-dir":outDir,"--json-conf":jsonconf}
        
        #merge with kwargs
        mergedOpts={**kwargs,**newOpts}
        
        status=self.runMikado("prepare",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if not status:
            print("Mikado prepare failed for:"+jsonconf)
            return ""
        
        #check if bam file exists
        if not checkPathsExists(outDir):
            return ""
        
        return outDir
        
        
        
    def runMikadoSerialise(self,jsonconf,blastTargets,orfs,xml,outDir="",verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper to run mikado serialise
        """
        #check input files exist
        if not check_files_exist(blastTargets,orfs,xml):
            print("Please check the input to mikado.")
            return ""
        if not outDir:
            outDir=os.getcwd()
        
        newOpts={"--json-conf":jsonconf,"--blast_targets":blastTargets,"--xml":xml,"--orfs":orfs,"--output-dir":outDir}
        
        #merge with kwargs
        mergedOpts={**kwargs,**newOpts}
        
        status=self.runMikado("serialise",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if not status:
            print("Mikado serialise failed for:"+jsonconf)
            return ""
        
        #check if bam file exists
        if not checkPathsExists(outDir):
            return ""
        
        return outDir
        
        
    def runMikadoPick(self,jsonconf,outDir="",verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper to run mikado pick
        """
        #check input files exist
        if not check_files_exist(jsonconf):
            print("Please check the input to mikado.")
            return ""
        if not outDir:
            outDir=os.getcwd()
        
        newOpts={"--json-conf":jsonconf}
        
        #merge with kwargs
        mergedOpts={**kwargs,**newOpts}
        
        status=self.runMikado("pick",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if not status:
            print("Mikado pick failed for:"+jsonconf)
            return ""
        
        #check if bam file exists
        if not checkPathsExists(outDir):
            return ""
        
        return outDir
        
        
    def runMikado(self,subCommand,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper to run mikado
        """
        
        #override existing arguments
        mergedArgsDict={**self.passedArgumentDict,**kwargs}
       
        mikado_Cmd=['mikado',subCommand]
        #add options
        mikado_Cmd.extend(parse_unix_args(self.validArgsList,mergedArgsDict))
                
        #print("Executing:"+" ".join(mergedArgsDict))
        
        #start ececution
        status=executeCommand(mikado_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            print_boldred("mikado failed")
        #return status
        return status
        
        
        
        
        
        
class Ribocode(RNASeqTools):
    def __init__(self,**kwargs):
        self.programName="RiboCode"
        self.depList=[self.programName]
        #check if program exists
        if not check_dependencies(self.depList):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        self.validArgsList=[]        
        self.passedArgumentDict=kwargs
        
        
    def runRibocode(self,gtf,genome,bam,l="no",outsuffix="ribocode_out",verbose=False,quiet=False,logs=True,objectid="NA"):
        """Wrapper to run ribocode in one step
        """
        
        #check input
        if not check_files_exist(gtf,genome,bam):
            print_boldred("Please check input files for Ribocode")
            return ""
        
        outDir=getFileDirectory(gtf)
        outFile=os.path.join(outDir,outsuffix)
        
        newOpts={"-g":gtf,"f":genome,"-r":bam,"-l":l,-o:outFile}
        
        ribocode_Cmd=['RiboCode_onestep']
        ribocode_Cmd.extend(parse_unix_args(self.validArgsList,newOpts))
        
        status=executeCommand(ribocode_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            print_boldred("ribocode failed")
            return ""
        
        return outFile
        
        
        

        
       
        
        
        
        
        
        
        
        
        
        
        