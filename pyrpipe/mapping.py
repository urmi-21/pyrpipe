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


class Kallisto(Aligner):
    """Kallisto constructor. Initialize kallisto parameters.
        """       
    def __init__(self,kallisto_index,**kwargs):
        super().__init__() 
        self.programName="kallisto"
        self.depList=[self.programName]        
        if not checkDep(self.depList):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        
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
        
        self.validArgsList=getListUnion(self.validArgsIndex,self.validArgsQuant,self.validArgsPseudo,self.validArgsh5dump)
        
        #initialize the passed arguments
        self.passedArgumentDict=kwargs
        
        #if index is passed, update the passed arguments
        if len(kallisto_index)>0 and checkFilesExists(kallisto_index):
            print("kallisto index is: "+kallisto_index)
            self.kallisto_index=kallisto_index
            self.passedArgumentDict['-i']=self.kallisto_index
        else:
            print("No kallisto index provided. Please use build_index() now to generate an index...")
            
    def build_kallisto_index(self,index_path,index_name,fasta,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """
        build kallisto index
        """
        
        #check input
        if not checkFilesExists(fasta):
            printBoldRed("{} does not exist. Exiting".format(fasta))
            return False
        
        #create out dir
        if not checkPathsExists(index_path):
            if not mkdir(indexPath):
                print("ERROR in building kallisto index. Failed to create index directory.")
                return False
            
        indexOut=os.path.join(index_path,index_name)
        newOpts={"--":(fasta,),"-i":indexOut}
        mergedOpts={**kwargs,**newOpts}
        
        #call salmon
        status=self.run_kallisto("index",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if status:
            #check if sam file is present in the location directory of sraOb
            if checkFilesExists(indexOut):
                self.kallisto_index=indexOut
                self.passedArgumentDict['-i']=self.kallisto_index
                printGreen("kallisto_index is:"+self.kallisto_index)
                return True
        else:
            printBoldRed("Failed to create kallisto index")
            return False
    
    def run_kallisto_quant(self,sraOb,outDir="",verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """
        run kallisto quant
        
        Returns
        -------
        string
            Path to salmon out directory
        """
        
        if not outDir:
            outDir=os.path.join(sraOb.location,"kallisto_out")
        
        
        
        if sraOb.layout == 'PAIRED':
            newOpts={"-o":outDir,"--":(sraOb.localfastq1Path,sraOb.localfastq2Path)}
        else:
            newOpts={"-o":outDir,"--single":"", "--":(sraOb.localfastqPath,)}
        
        
        #add input files to kwargs, overwrite kwargs with newOpts
        mergedOpts={**kwargs,**newOpts}
        
        #call salmon
        status=self.run_kallisto("quant",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if status:
            #check if sam file is present in the location directory of sraOb
            if checkFilesExists(os.path.join(outDir,"abundance.tsv")):
                return outDir
        
        printBoldRed("kallisto quant failed")
        return ""
        
    
    
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
    
    def checkIndex(self):
        return checkFilesExists(self.kallisto_index)
            


class Salmon(Aligner):
    """Salmon constructor. Initialize kallisto parameters.
    """       
    def __init__(self,salmon_index,**kwargs):    
        super().__init__() 
        self.programName="salmon"
        self.depList=[self.programName]        
        if not checkDep(self.depList):
            raise Exception("ERROR: "+ self.programName+" not found.")
        
        
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

        self.validArgsList=getListUnion(self.validArgsIndex,self.validArgsQuantReads,self.validArgsQuantAlign,self.validArgsQuantMerge)
        
        #initialize the passed arguments
        self.passedArgumentDict=kwargs
        
        #if index is passed, update the passed arguments
        if len(salmon_index)>0 and check_salmon_index(salmon_index):
            print("salmon index is: "+salmon_index)
            self.salmon_index=salmon_index
            self.passedArgumentDict['-i']=self.salmon_index
        else:
            print("No salmon index provided. Please build index now to generate an index...")
            
            
            
    def build_salmon_index(self,index_path,index_name,fasta,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """
        build salmon index
        """
        
        #check input
        if not checkFilesExists(fasta):
            printBoldRed("{} does not exist. Exiting".format(fasta))
            return False
        #create out dir
        if not checkPathsExists(index_path):
            if not mkdir(indexPath):
                print("ERROR in building hisat2 index. Failed to create index directory.")
                return False
        indexOut=os.path.join(index_path,index_name)
        newOpts={"-t":fasta,"-i":indexOut}
        mergedOpts={**kwargs,**newOpts}
        
        #call salmon
        status=self.run_salmon("index",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if status:
            #check if sam file is present in the location directory of sraOb
            #if checkFilesExists(os.path.join(indexOut,"versionInfo.json")): #not sure if this is reliable
            if checkPathsExists(indexOut):
                self.salmon_index=indexOut
                self.passedArgumentDict['-i']=self.salmon_index
                printGreen("salmon index is:"+self.salmon_index)
                return True
        
        printBoldRed("Failed to create salmon index")
        return False
        
        
    
    def run_salmon_quant(self,sraOb,outDir="",libType="A",verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """
        run salmon quant
        
        Returns
        -------
        string
            Path to salmon out directory
        """
        
        if not outDir:
            outDir=os.path.join(sraOb.location,"salmon_out")
        
        
        
        if sraOb.layout == 'PAIRED':
            newOpts={"-o":outDir,"-l":libType,"-1":sraOb.localfastq1Path,"-2":sraOb.localfastq2Path}
        else:
            newOpts={"-o":outDir,"-l":libType,"-r":sraOb.localfastqPath}
        
        
        #add input files to kwargs, overwrite kwargs with newOpts
        mergedOpts={**kwargs,**newOpts}
        
        #call salmon
        status=self.run_salmon("quant",verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if status:
            #check if sam file is present in the location directory of sraOb
            if checkFilesExists(os.path.join(outDir,"quant.sf")):
                return outDir
        
        printBoldRed("salmon quant failed")
        return ""
        
        
        
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
                raise Exception("ERROR: Invalid salmon index. Please run build index to generate an index.")
        
        #override existing arguments
        mergedArgsDict={**self.passedArgumentDict,**kwargs}
            
        salmon_Cmd=['salmon',subcommand]
        salmon_Cmd.extend(parseUnixStyleArgs(self.validArgsList,mergedArgsDict))
        
        #start ececution
        status=executeCommand(salmon_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,command_name=" ".join(salmon_Cmd[0:2]))
        if not status:
            printBoldRed("salmon failed")
        return status 

    def checkIndex(self):
        if hasattr(self,'salmon_index'):
            return checkSalmonIndex(self.salmon_index)
        return False


         