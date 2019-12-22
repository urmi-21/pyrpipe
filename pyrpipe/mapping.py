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

class Hisat2(Aligner):
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
            
        
        
            
    def buildHisat2Index(self,indexPath,indexName,*args,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
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
        #print("Executing:"+str(" ".join(hisat2Build_Cmd)))
        
        #start ececution
        status=executeCommand(hisat2Build_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
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
        
        
    def performAlignment(self,sraOb,outSamSuffix="_hisat2",verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
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
        status=self.runHisat2(verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if status:
            #check if sam file is present in the location directory of sraOb
            if checkFilesExists(outSamFile):
                return outSamFile
        else:
            return ""
            
        
    def runHisat2(self,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
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
        cmdStatus=executeCommand(hisat2_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not cmdStatus:
            print("hisat2 failed:"+" ".join(hisat2_Cmd))
     
        #return status
        return cmdStatus
        
        
    
    def checkHisat2Index(self):
        if hasattr(self,'hisat2Index'):
            return(checkHisatIndex(self.hisat2Index))
        else:
            return False

    



class Star(Aligner):
    def __init__(self,starIndex="",**kwargs):
        """STAR constructor. Initialize star's index and other parameters.
        """
        super().__init__() 
        self.programName="STAR"
        
        self.depList=[self.programName]        
        #check if hisat2 exists
        if not checkDep(self.depList):
            raise Exception("ERROR: "+ self.programName+" not found.")
            
        self.validArgsList=['--help','--parametersFiles','--sysShell','--runMode','--runThreadN','--runDirPerm','--runRNGseed','--quantMode','--quantTranscriptomeBAMcompression','--quantTranscriptomeBan','--twopassMode','--twopass1readsN',
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
        if len(starIndex)>0 and checkStarIndex(starIndex):
            print("STAR index is: "+starIndex)
            self.starIndex=starIndex
            self.passedArgumentDict['--genomeDir']=self.starIndex
        else:
            print("No STAR index provided. Please build index now to generate an index using buildStarIndex()....")
            
    
    #STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ./starIndex --genomeFastaFiles /home/usingh/work/urmi/hoap/test/hisatYeast/S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113.fsa
    def buildStarIndex(self,indexPath,*args,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Build a star index with given parameters and saves the new index to self.starIndex.
        Parameters
        ----------
        arg1: string
            Path where the index will be created
        
        arg2: tuple
            Path to reference input files
        arg3: dict
            Parameters for the star command
        
        Returns
        -------
        bool:
            Returns status of star command
        """
        if len(args)<1:
            printBoldRed("Please provide input fasta file to build STAR index")
            return ""
        
        print("Building STAR index...")
        
        #create path if doesnt exists
        if not checkPathsExists(indexPath):
            if not mkdir(indexPath):
                raise Exception("Error creating STAR index. Exiting.")
                return False
        
        #add runMode
        newOpts={"--runMode":"genomeGenerate","--genomeDir":indexPath,"--genomeFastaFiles":" ".join(args)}
        
        mergedOpts={**kwargs,**newOpts}
        
        starbuild_Cmd=['STAR']
        starbuild_Cmd.extend(parseUnixStyleArgs(self.validArgsList,mergedOpts))
        
        #execute command
        status=executeCommand(starbuild_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        
        
        if status:
            print("Star build finished")
            #check if sam file is present in the location directory of sraOb
            if checkPathsExists(indexPath):
                #update object's index
                self.starIndex=indexPath
                self.passedArgumentDict['--genomeDir']=self.starIndex
                if self.checkstarIndex():
                    return True
        else:
            return False
        
 
            
    def performAlignment(self,sraOb,outSamSuffix="_star",verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Function to perform alignment using self object and the provided sraOb.
        All star output will be written to the sraOb directory by default.
        
        Parameters
        ----------
        arg1: SRA object
            An object of type SRA. The path to fastq files will be obtained from this object.
        arg2: string
            Suffix for the output file
        arg3: dict
            Options to pass to hisat2.
            
        Returns
        -------
        string:
            path to the output dir
        """
        outDir=sraOb.location
        
        #find layout and fq file paths
        if sraOb.layout == 'PAIRED':
            newOpts={"--readFilesIn":sraOb.localfastq1Path+" "+sraOb.localfastq2Path}
        else:
            newOpts={"--readFilesIn":sraOb.localfastqPath}
        
        #add out dir
        newOpts["--outFileNamePrefix"]=outDir+"/"
        
        #add index
        
        
        #add input files to kwargs, overwrite kwargs with newOpts
        mergedOpts={**kwargs,**newOpts}
        
        #call star
        status=self.runStar(verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
                
        
        if status:
            print("Star finished")
            #check if sam file is present in the location directory of sraOb
            if checkPathsExists(outDir):
                
                return outDir
        else:
            return ""
        
        
    
    def runStar(self,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper for running star.
        The self.starIndex index used.
        
        Parameters
        ----------
        arg1: dict
            arguments to pass to star. This will override parametrs already existing in the self.passedArgumentList list but NOT replace all of them.
            
        Returns
        -------
        bool:
                Returns the status of star. True is passed, False if failed.
        """
        
        #check for a valid index
        if not self.checkstarIndex():
            raise Exception("ERROR: Invalid star index. Please run build index to generate an index.")
            
        #override existing arguments
        mergedArgsDict={**self.passedArgumentDict,**kwargs}
       
        star_Cmd=['STAR']
        #add options
        star_Cmd.extend(parseUnixStyleArgs(self.validArgsList,mergedArgsDict))        
        
        #execute command
        cmdStatus=executeCommand(star_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        
        if not cmdStatus:
            print("STAR failed:"+" ".join(star_Cmd))
     
        #return status
        return cmdStatus
    
    
    def checkstarIndex(self):
        if hasattr(self,'starIndex'):
            return(checkStarIndex(self.starIndex))
        else:
            return False
            



class Bowtie2(Aligner):
    def __init__(self,bowtie2Index,**kwargs):
        """Bowtie2 constructor. Initialize bowtie2 index and other parameters.
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
        
        
    
    
    def performAlignment(self,sraOb,outSamSuffix="_bt2",overwrite=True,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
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
        
        status=self.runBowTie2(verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if status:
            #check if sam file is present in the location directory of sraOb
            if checkFilesExists(outFile):
                return outFile
        else:
            return ""
        
        
        
    
    def runBowTie2(self,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
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
        
        #print("Executing:"+" ".join(bowtie2_Cmd))
        
        #start ececution
        status=executeCommand(bowtie2_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            printBoldRed("bowtie2 failed")
        return status
    
    
    def checkIndex(self):
        if hasattr(self,'bowtie2Index'):
            return(checkBowtie2Index(self.bowtie2Index))
        else:
            return False

#TODO
class Kallisto(Aligner):
    """Kallisto constructor. Initialize kallisto parameters.
        """       
    def __init__(self,kallisto_index,**kwargs):
        super().__init__() 
        self.programName="kallisto"
        self.depList=[self.programName]        
        if not checkDep(self.depList):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        self.validArgsList=[]
        
        #initialize the passed arguments
        self.passedArgumentDict=kwargs
        
        #if index is passed, update the passed arguments
        if len(kallisto_index)>0 and checkKallistoIndex(kallisto_index):
            print("kallisto index is: "+kallisto_index)
            self.kallisto_index=kallisto_index
            self.passedArgumentDict['-i']=self.kallisto_index
        else:
            print("No Bowtie2 index provided. Please build index now to generate an index...")
            
    def build_kallisto_index(self,index_path,index_name,*args,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """
        build kallisto index
        """
        pass
    
    def run_kallisto_quant(self,index_path,index_name,*args,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """
        run kallisto quant
        """
        pass
    def run_kallisto(self,subcommand,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper for running kallisto.
        
        ----------
        arg1: dict
            arguments to pass to bowtie2. This will override parametrs already existing in the self.passedArgumentList list but NOT replace them.
            
        Returns
        -------
        bool:
                Returns the status of bowtie2. True is passed, False if failed.
        """
        
        #check for a valid index
        if subcommand!="index":
            if not self.checkIndex():
                raise Exception("ERROR: Invalid kallisto index. Please run build index to generate an index.")
        
        #override existing arguments
        mergedArgsDict={**self.passedArgumentDict,**kwargs}
            
        kallisto_Cmd=['kallisto',subcommand]
        kallisto_Cmd.extend(parseUnixStyleArgs(self.validArgsList,mergedArgsDict))
        
        #start ececution
        status=executeCommand(kallisto_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,command_name=" ".join(kallisto_Cmd[0:2]))
        if not status:
            printBoldRed("kallisto failed")
        return status       
            
            
class Salmon(Aligner):
    """Salmon constructor. Initialize kallisto parameters.
        """       
    def __init__(self,salmon_index,**kwargs):    
        super().__init__() 
        self.programName="salmon"
        self.depList=[self.programName]        
        if not checkDep(self.depList):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        self.validArgsList=[]
        
        #initialize the passed arguments
        self.passedArgumentDict=kwargs
        
        #if index is passed, update the passed arguments
        if len(salmon_index)>0 and check_salmon_index(salmon_index):
            print("salmon index is: "+salmon_index)
            self.salmon_index=salmon_index
            self.passedArgumentDict['-i']=self.salmon_index
        else:
            print("No salmon index provided. Please build index now to generate an index...")
            
            
            
    def build_salmon_index(self,index_path,index_name,*args,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """
        build salmon index
        """
        pass
    
    def run_salmon_quant(self,index_path,index_name,*args,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """
        run salmon quant
        """
        pass
    def run_salmon(self,subcommand,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper for running kallisto.
        
        ----------
        arg1: dict
            arguments to pass to bowtie2. This will override parametrs already existing in the self.passedArgumentList list but NOT replace them.
            
        Returns
        -------
        bool:
                Returns the status of bowtie2. True is passed, False if failed.
        """
        
        #check for a valid index
        if subcommand!="index":
            if not self.checkIndex():
                raise Exception("ERROR: Invalid kallisto index. Please run build index to generate an index.")
        
        #override existing arguments
        mergedArgsDict={**self.passedArgumentDict,**kwargs}
            
        salmon_Cmd=['salmon',subcommand]
        salmon_Cmd.extend(parseUnixStyleArgs(self.validArgsList,mergedArgsDict))
        
        #start ececution
        status=executeCommand(salmon_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,command_name=" ".join(salmon_Cmd[0:2]))
        if not status:
            printBoldRed("salmon failed")
        return status 




         