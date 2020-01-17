#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 19:53:42 2019

@author: usingh
contains classes of RNA-Seq alignment programs
"""

from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe
import os

class Aligner:
    """This is an abstract class for alignment programs.
    """
    def __init__(self,index=""):
        self.category="Aligner"
        self.passedArgumentDict={}
        self.index=index
        
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
    """This class represents hisat2 program.
    
       Parameters
       ----------       
       
       hisat2_index: string
            path to the histat2 index. This index is stored with the object and will be used when hisat is invoked using this object.
            
       **kwargs: dict
            parameters passed to the hisat2 program. These parameters could be overridden later when running hisat.
            
    Attributes
    ----------
    
    """ 
    def __init__(self,hisat2_index="",**kwargs):
        
        super().__init__() 
        self.programName="hisat2"
        #check if hisat2 exists
        if not pe.check_dependencies([self.programName]):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        self.valid_args=['-x','-1','-2','-U','--sra-acc','-S','-q','--qseq','-f','-r','-c','-s',
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
        if len(hisat2_index)>0 and pu.check_hisatindex(hisat2_index):
            print("HISAT2 index is: "+hisat2_index)
            self.hisat2_index=hisat2_index
            self.passedArgumentDict['-x']=self.hisat2_index
            self.index=self.hisat2_index
        else:
            print("No Hisat2 index provided. Please build index now to generate an index using build_Index()....")
            
        
        
            
    def build_index(self,index_path,index_name,*args,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Build a hisat index with given parameters and saves the new index to self.hisat2_index.
        
        Parameters
        ----------
        
        index_path: string
            Path where the index will be created
            
        index_name: string
            A name for the index
            
        args: tuple
            Path to reference input files
            
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
        
        #check input references
        if len(args)<1:
            pu.print_boldred("No reference sequence provided to hisat2-build. Exiting")
            return False
        
        if not pu.check_files_exist(*args):
            pu.print_boldred("Please check input reference sequences provided to hisat2-build. Exiting")
            return False
            
        overwrite=True
        print("Building hisat index...")
        
        hisat2Buildvalid_args=['-c','--large-index','-a','-p','--bmax','--bmaxdivn','--dcv','--nodc','-r','-3','-o',
                                  '-t','--localoffrate','--localftabchars','--snp','--haplotype','--ss','--exon',
                                  '--seed','-q','-h','--usage','--version']
        #create the out dir
        if not pu.check_paths_exist(index_path):
            if not pu.mkdir(index_path):
                print("ERROR in building hisat2 index. Failed to create index directory.")
                return False
        
        if not overwrite:
            #check if files exists
            if pu.check_hisatindex(os.path.join(index_path,index_name)):
                print("Hisat2 index with same name already exists. Exiting...")
                return True
        
        hisat2Build_Cmd=['hisat2-build']
        #add options
        hisat2Build_Cmd.extend(pu.parse_unix_args(hisat2Buildvalid_args,kwargs))
        #add input files
        hisat2Build_Cmd.append(str(",".join(args)))
        #add dir/basenae
        hisat2Build_Cmd.append(os.path.join(index_path,index_name))
        #print("Executing:"+str(" ".join(hisat2Build_Cmd)))
        
        #start ececution
        status=pe.execute_command(hisat2Build_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            pu.print_boldred("hisatBuild failed")
            return False
        
        #check index files
        if not pu.check_hisatindex(os.path.join(index_path,index_name)):
            pu.print_boldred("hisatBuild failed")
            return False
        
        #set the index path
        self.hisat2_index=os.path.join(index_path,index_name)
        self.passedArgumentDict['-x']=self.hisat2_index
        
        #return status
        return True
        
        
    def perform_alignment(self,sra_object,out_suffix="_hisat2",verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Function to perform alignment using sra_object.
        
        Parameters
        ----------
        
        sra_object SRA object
            An object of type SRA. The path to fastq files will be obtained from this object.
        out_suffix: string
            Suffix for the output sam file
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to hisat2. This will override the existing options in self.passed_args_dict (only replace existing arguments and not replace all the arguments).
        """
        
        
        #create path to output sam file
        outSamFile=os.path.join(sra_object.location,sra_object.srr_accession+out_suffix+".sam")
        
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
        if sra_object.layout == 'PAIRED':
            newOpts={"-1":sra_object.localfastq1Path,"-2":sra_object.localfastq2Path,"-S":outSamFile}
        else:
            newOpts={"-U":sra_object.localfastqPath,"-S":outSamFile}
        
        #add input files to kwargs, overwrite kwargs with newOpts
        mergedOpts={**kwargs,**newOpts}
        
        #call run_hisat2
        status=self.run_hisat2(verbose=verbose,quiet=quiet,logs=logs,objectid=sra_object.srr_accession,**mergedOpts)
        
        if status:
            #check if sam file is present in the location directory of sra_object
            if pu.check_files_exist(outSamFile):
                return outSamFile
        else:
            return ""
            
        
    def run_hisat2(self,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper for running hisat2.
        
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

        :return: Returns the status of hisat2. True is passed, False if failed.
        :rtype: bool
        """
        
        #check for a valid index
        if not self.check_index():
            raise Exception("ERROR: Invalid HISAT2 index. Please run build index to generate an index.")
            
        #override existing arguments
        mergedArgsDict={**self.passedArgumentDict,**kwargs}
       
        hisat2_Cmd=['hisat2']
        #add options
        hisat2_Cmd.extend(pu.parse_unix_args(self.valid_args,mergedArgsDict))        
        
        #execute command
        cmd_status=pe.execute_command(hisat2_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not cmd_status:
            print("hisat2 failed:"+" ".join(hisat2_Cmd))
     
        #return status
        return cmd_status
        
        
    
    def check_index(self):
        if hasattr(self,'hisat2_index'):
            return(pu.check_hisatindex(self.hisat2_index))
        else:
            return False

    



class Star(Aligner):
    """This class represents STAR program.
    
       Parameters
       ----------
       
       star_index: string
            path to a star index. This index will be used when star is invoked using this object.
       kwargs: dict
            parameters passed to the star program. These parameters could be overridden later when running star.
            
    Attributes
    ----------
    """ 
    def __init__(self,star_index="",**kwargs):
        
        super().__init__() 
        self.programName="STAR"
        
        self.dep_list=[self.programName]        
        #check if star exists
        if not pe.check_dependencies(self.dep_list):
            raise Exception("ERROR: "+ self.programName+" not found.")
            
        self.valid_args=['--help','--parametersFiles','--sysShell','--runMode','--runThreadN','--runDirPerm','--runRNGseed','--quantMode','--quantTranscriptomeBAMcompression','--quantTranscriptomeBan','--twopassMode','--twopass1readsN',
                            '--genomeDir','--genomeLoad','--genomeFastaFiles','--genomeChrBinNbits','--genomeSAindexNbases','--genomeSAsparseD','--genomeSuffixLengthMax','--genomeChainFiles','--genomeFileSizes',
                            '--sjdbFileChrStartEnd','--sjdbGTFfile','--sjdbGTFchrPrefix','--sjdbGTFfeatureExon','--sjdbGTFtagExonParentTranscript','--sjdbGTFtagExonParentGene','--sjdbOverhang','--sjdbScore','--sjdbInsertSave',
                            '--inputBAMfile','--readFilesIn','--readFilesCommand','--readMapNumber','--readMatesLengthsIn','--readNameSeparator','--clip3pNbases','--clip5pNbases','--clip3pAdapterSeq','--clip3pAdapterMMp','--clip3pAfterAdapterNbases',
                            '--limitGenomeGenerateRAM','--limitIObufferSize','--limitOutSAMoneReadBytes','--limitOutSJoneRead','--limitOutSJcollapsed','--limitBAMsortRAM ','--limitSjdbInsertNsj','--outFileNamePrefix','--outTmpDir','--outTmpKeep',
                            '--outStd','--outReadsUnmapped','--outQSconversionAdd','--outMultimapperOrder','--outSAMtype','--outSAMmode','--outSAMstrandField','--outSAMattributes','--outSAMattrIHstart','--outSAMunmapped','--outSAMorder',
                            '--outSAMprimaryFlag','--outSAMreadID','--outSAMmapqUnique','--outSAMflagOR','--outSAMflagAND','--outSAMattrRGline','--outSAMheaderHD','--outSAMheaderPG','--outSAMheaderCommentFile','--outSAMfilter','--outSAMmultNmax',
                            '--outBAMcompression','--outBAMsortingThreadN','--bamRemoveDuplicatesType','--bamRemoveDuplicatesMate2basesN','--outWigType','--outWigStrand','--outWigReferencesPrefix','--outWigNorm','--outFilterType',
                            '--outFilterMultimapScoreRange','--outFilterMultimapNmax','--outFilterMismatchNmax','--outFilterMismatchNoverLmax','--outFilterMismatchNoverReadLmax','--outFilterScoreMin','--outFilterScoreMinOverLread',
                            '--outFilterMatchNmin','--outFilterMatchNminOverLread','--outFilterIntronMotifs','--outSJfilterReads','--outSJfilterOverhangMin','--outSJfilterCountUniqueMin','--outSJfilterCountTotalMin','--outSJfilterDistToOtherSJmin',
                            '--outSJfilterIntronMaxVsReadN','--scoreGap','--scoreGapNoncan','--scoreGapGCAG ','--scoreGapATAC','--scoreGenomicLengthLog2scale','--scoreDelOpen','--scoreDelBase','--scoreInsOpen','--scoreInsBase','--scoreStitchSJshift',
                            '--seedSearchStartLmax','--seedSearchStartLmaxOverLread','--seedSearchLmax','--seedMultimapNmax','--seedPerReadNmax','--seedPerWindowNmax','--seedNoneLociPerWindow','--alignIntronMin','--alignIntronMax','--alignMatesGapMax',
                            '--alignSJoverhangMin','--alignSJstitchMismatchNmax','--alignSJDBoverhangMin','--alignSplicedMateMapLmin','--alignSplicedMateMapLminOverLmate','--alignWindowsPerReadNmax','--alignTranscriptsPerWindowNmax','--alignTranscriptsPerReadNmax',
                            '--alignEndsType','--alignEndsProtrude','--alignSoftClipAtReferenceEnds','--winAnchorMultimapNmax','--winBinNbits','--winAnchorDistNbins','--winFlankNbins','--winReadCoverageRelativeMin','--winReadCoverageBasesMin',
                            '--chimOutType','--chimSegmentMin','--chimScoreMin','--chimScoreDropMax','--chimScoreSeparation','--chimScoreJunctionNonGTAG','--chimJunctionOverhangMin','--chimSegmentReadGapMax','--chimFilter','--chimMainSegmentMultNmax']
                
       
        #initialize the passed arguments
        self.passedArgumentDict=kwargs
        
        #if index is passed, update the passed arguments
        if len(star_index)>0 and pu.check_starindex(star_index):
            print("STAR index is: "+star_index)
            self.star_index=star_index
            self.passedArgumentDict['--genomeDir']=self.star_index
        else:
            print("No STAR index provided. Please build index now to generate an index using build_index()....")
            
    
    
    def build_index(self,index_path,*args,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Build a star index with given parameters and saves the new index to self.star_index.
        
        Parameters
        ----------
        
        index_path: string
            Path where the index will be created
        
        args: tuple
            Path to reference input files
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        
        kwargs: dict
            Parameters for the star command

        :return: Returns status of star command
        :rtype: bool
        """
        
        #if index already exists then exit
        if pu.check_starindex(index_path):
            pu.print_green("STAR index already exists. Using it...")
            self.star_index=index_path
            self.passedArgumentDict['--genomeDir']=self.star_index
            return True
            
        
        #check input files
        if len(args)<1:
            pu.print_boldred("Please provide input fasta file to build STAR index")
            return False
        
        if not pu.check_files_exist(*args):
            raise Exception("Please check input to star index")
            return False
        
        #create path if doesnt exist
        if not pu.check_paths_exist(index_path):
            if not pu.mkdir(index_path):
                raise Exception("Error creating STAR index. Exiting.")
                return False
        
        #add runMode
        newOpts={"--runMode":"genomeGenerate","--genomeDir":index_path,"--genomeFastaFiles":" ".join(args)}
        
        mergedOpts={**kwargs,**newOpts}
        
        starbuild_Cmd=['STAR']
        starbuild_Cmd.extend(pu.parse_unix_args(self.valid_args,mergedOpts))
        
        #execute command
        status=pe.execute_command(starbuild_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        
        
        if status:
            if pu.check_paths_exist(index_path):
                #update object's index
                self.star_index=index_path
                self.passedArgumentDict['--genomeDir']=self.star_index
                if self.check_index():
                    return True
        else:
            return False
        
 
            
    def perform_alignment(self,sra_object,out_suffix="_star",out_dir="",verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Function to perform alignment using star and the provided SRA object.
        All star output will be written to the sra_object directory by default.
        
        Parameters
        ----------
        
        sra_object: SRA object
            An object of type SRA. The path to fastq files will be obtained from this object.
        out_suffix: string
            Suffix for the output file
        out_dir: str
            outout directory default: sra_object.location
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        kwargs: dict
            Options to pass to STAR. This will override the existing options in self.passed_args_dict (only replace existing arguments and not replace all the arguments).
        

        :return: Return the path to the output dir
        :rtype: string
        """
        
        if not out_dir:
            out_dir=sra_object.location
        else:
            #create out_dir if not exists
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
        
        #find layout and fq file paths
        if sra_object.layout == 'PAIRED':
            newOpts={"--readFilesIn":sra_object.localfastq1Path+" "+sra_object.localfastq2Path}
        else:
            newOpts={"--readFilesIn":sra_object.localfastqPath}
        
        #add out dir
        newOpts["--outFileNamePrefix"]=out_dir+"/"
        
               
        #add input files to kwargs, overwrite kwargs with newOpts
        mergedOpts={**kwargs,**newOpts}
        
        #call star
        status=self.run_star(verbose=verbose,quiet=quiet,logs=logs,objectid=sra_object.srr_accession,**mergedOpts)
                
        
        if status:
            print("Star finished")
            if pu.check_paths_exist(out_dir):
                return out_dir
        else:
            return ""
        
    
    def run_star(self,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper for running star.
        The self.star_index index used.
        
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
            Options to pass to stringtie. This will override the existing options in self.passed_args_dict (only replace existing arguments and not replace all the arguments).
        kwargs: dict
            arguments to pass to star. This will override parametrs already existing in the self.passedArgumentList list but NOT replace all of them.
            
        :return: Returns the status of star. True is passed, False if failed.
        :rtype: bool
        """
        
        #check for a valid index
        if not self.check_index():
            raise Exception("ERROR: Invalid star index. Please run build index to generate an index.")
            
        #override existing arguments
        mergedArgsDict={**self.passedArgumentDict,**kwargs}
       
        star_cmd=['STAR']
        #add options
        star_cmd.extend(pu.parse_unix_args(self.valid_args,mergedArgsDict))        
        
        #execute command
        cmd_status=pe.execute_command(star_cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        
        if not cmd_status:
            print("STAR failed:"+" ".join(star_cmd))
     
        #return status
        return cmd_status
    
    
    def check_index(self):
        if hasattr(self,'star_index'):
            return(pu.check_starindex(self.star_index))
        else:
            return False
            



class Bowtie2(Aligner):
    """This class represents bowtie2 program.
    
       Parameters
       ----------
       
       bowtie2_index: string
            path to a bowtie2 index. This index will be used when bowtie2 is invoked using this object.
       **kwargs: dict
            parameters passed to the bowtie2 program. These parameters could be overridden later when running bowtie2.
            
       Attributes
       ----------
    """ 
    def __init__(self,bowtie2_index,**kwargs):
        """Bowtie2 constructor. Initialize bowtie2 index and other parameters.
        """       
        
        super().__init__() 
        self.programName="bowtie2"
        self.dep_list=[self.programName]        
        if not pe.check_dependencies(self.dep_list):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        self.valid_args=['-x','-1','-2','-U','--interleaved','-S','-b','-q','--tab5','--tab6','--qseq','-f','-r','-F','-c','-s','-u','-5','-3',
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
        if len(bowtie2_index)>0 and pu.check_bowtie2index(bowtie2_index):
            print("Bowtie2 index is: "+bowtie2_index)
            self.bowtie2_index=bowtie2_index
            self.passedArgumentDict['-x']=self.bowtie2_index
        else:
            print("No Bowtie2 index provided. Please build index now to generate an index...")
        
        
    def build_index(self,index_path,index_name,*args,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Build a bowtie2 index with given parameters and saves the new index to self.bowtie2_index.
        
        Parameters
        ----------
        
        index_path: string
            Path where the index will be created
        index_name: string
            A name for the index
        arg3: tuple
            Path to reference input files
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
        if len(args)<1:
            pu.print_boldred("No reference sequence provided to bowtie2-build. Exiting")
            return False
        
        if not pu.check_files_exist(*args):
            pu.print_boldred("Please check input reference sequences provided to bowtie2-build. Exiting")
            return False
            
        overwrite=True
        
        
        bowtie2_build_args=['-f','-c','--large-index','--debug','--sanitized','--verbose','-a',
                            '--noauto','-p','--packed','--bmax','--bmaxdivn','--dcv','--nodc',
                            '-r','--noref','-3','--justref','-o','--offrate','-t','--ftabchars',
                            '--threads','--seed','-q','--quiet','-h','--help','--usage','--version']
        #create the out dir
        if not pu.check_paths_exist(index_path):
            if not pu.mkdir(index_path):
                print("ERROR in building bowtie2 index. Failed to create index directory.")
                return False
        
        if not overwrite:
            #check if files exists
            if pu.check_bowtie2index(os.path.join(index_path,index_name)):
                print("bowtie2 index with same name already exists. Exiting...")
                return True
            
        bowtie2Build_Cmd=['bowtie2-build']
        #add options
        bowtie2Build_Cmd.extend(pu.parse_unix_args(bowtie2_build_args,kwargs))
        #add input files
        bowtie2Build_Cmd.append(str(",".join(args)))
        #add dir/basenae
        bowtie2Build_Cmd.append(os.path.join(index_path,index_name))
        #print("Executing:"+str(" ".join(hisat2Build_Cmd)))
        
        #start ececution
        status=pe.execute_command(bowtie2Build_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            pu.print_boldred("bowtie2-build failed")
            return False
        
        #check index files
        if not pu.check_bowtie2index(os.path.join(index_path,index_name)):
            pu.print_boldred("bowtie2-build failed")
            return False
        
        #set the index path
        self.bowtie2_index=os.path.join(index_path,index_name)
        self.passedArgumentDict['-x']=self.bowtie2_index
        
        #return status
        return True
        
    
    def perform_alignment(self,sra_object,out_suffix="_bt2",out_dir="",overwrite=True,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Function to perform alignment using self object and the provided sra_object.
        
        Parameters
        ----------
        
        sra_object: SRA object
            An object of type SRA. The path to fastq files will be obtained from this object.
        out_suffix: string
            Suffix for the output sam file
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
            out_dir=sra_object.location
        else:
            if not pu.check_paths_exist(out_dir):
                pu.mkdir(out_dir)
                
        #create path to output sam file
        outFile=os.path.join(out_dir,sra_object.srr_accession+out_suffix+".sam")
                    
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
        if sra_object.layout == 'PAIRED':
            newOpts={"-1":sra_object.localfastq1Path,"-2":sra_object.localfastq2Path,"-S":outFile}
        else:
            newOpts={"-U":sra_object.localfastqPath,"-S":outFile}
        
        #add input files to kwargs, overwrite kwargs with newOpts
        mergedOpts={**kwargs,**newOpts}
        
        status=self.run_bowtie2(verbose=verbose,quiet=quiet,logs=logs,objectid=sra_object.srr_accession,**mergedOpts)
        
        if status:
            #check if sam file is present in the location directory of sra_object
            if pu.check_files_exist(outFile):
                return outFile
        else:
            return ""
        
        
        
    
    def run_bowtie2(self,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper for running bowtie2.
        
        
        
        verbose: bool
            Print stdout and std error
        quiet: bool
            Print nothing
        logs: bool
            Log this command to pyrpipe logs
        objectid: str
            Provide an id to attach with this command e.g. the SRR accession. This is useful for debugging, benchmarking and reports.
        
        kwargs: dict
            arguments to pass to bowtie2. This will override parametrs already existing in the self.passedArgumentList list but NOT replace them.
            
        :return: Returns the status of bowtie2. True is passed, False if failed.
        :rtype: bool
        """
        
        #check for a valid index
        if not self.check_index():
            raise Exception("ERROR: Invalid Bowtie2 index. Please run build index to generate an index.")
        
        #override existing arguments
        mergedArgsDict={**self.passedArgumentDict,**kwargs}
            
        bowtie2_cmd=['bowtie2']
        bowtie2_cmd.extend(pu.parse_unix_args(self.valid_args,mergedArgsDict))
        
        #print("Executing:"+" ".join(bowtie2_cmd))
        
        #start ececution
        status=pe.execute_command(bowtie2_cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            pu.print_boldred("bowtie2 failed")
        return status
    
    
    def check_index(self):
        """Function to check bowtie index.
        Returns True is index exist on disk.
        """
        if hasattr(self,'bowtie2_index'):
            return(pu.check_bowtie2index(self.bowtie2_index))
        return False



         
