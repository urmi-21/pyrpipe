#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:21:01 2019

@author: usingh
"""

from pyrpipe.pyrpipe_utils import *
from pyrpipe.pyrpipe_engine import *

class Assembly:
    def __init__(self):
        self.category="Assembler"
    def performAssembly():
        pass

class Stringtie(Assembly):
    def __init__(self,referenceGTF="",**kwargs):
        """Stringtie constructor. Initialize stringtie parameters.
        
        Parameters
        ----------
        arg1: string
            Path to the reference gtf file.
        arg2: dict
            Options passed to stringtie. These could be overridden later when executing stringtie.
        """
        super().__init__()
        self.programName="stringtie"
        #check if stringtie exists
        if not checkDep([self.programName]):
            raise Exception("ERROR: "+ self.programName+" not found.")
        self.validArgsList=['-G','--version','--conservative','--rf','--fr','-o','-l',
                            '-f','-L','-m','-a','-j','-t','-c','-s','-v','-g','-M',
                            '-p','-A','-B','-b','-e','-x','-u','-h','--merge','-F','-T','-i']
        
        #keep the passed arguments
        self.passedArgumentDict=kwargs
        
        #check the reference GTF
        if len(referenceGTF)>0 and checkFilesExists(referenceGTF):
            self.referenceGTF=referenceGTF
            self.passedArgumentDict['-G']=referenceGTF
        
    def performAssembly(self,inputBAM,outFileSuffix="_stringtie",overwrite=True,**kwargs):
        """Function to run stringtie using BAM file.
                
        Parameters
        ----------
        arg1: string
            path to bam file
        arg2: string
            Suffix for the output gtf file
        arg3: bool
            Overwrite if output file already exists.
        arg4: dict
            Options to pass to stringtie. This will override the existing options self.passedArgumentDict (only replace existing arguments and not replace all the arguments).
            
        Returns
        -------
        string
            path to output GTF file
        """
        
        #create path to output file
        fname=getFileBaseName(inputBAM)
        outDir=getFileDirectory(inputBAM)
        outGtfFile=os.path.join(outDir,fname+outFileSuffix+".gtf")
        
        """
        Handle overwrite
        """
        if not overwrite:
            #check if file exists. return if yes
            if os.path.isfile(outGtfFile):
                print("The file "+outGtfFile+" already exists. Exiting..")
                return outGtfFile
        
        #Add output file name and input bam
        newOpts={"-o":outGtfFile,"--":(inputBAM,)}
        mergedOpts={**kwargs,**newOpts}
        
        #call stringtie
        status=self.runStringtie(**mergedOpts)
        
        if status:
            #check if sam file is present in the location directory of sraOb
            if checkFilesExists(outGtfFile):
                return outGtfFile
        else:
            return ""
        
    def performStringtieMerge(self,*args,outFileSuffix="_stringtieMerge",overwrite=True,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Function to run stringtie merge.
        Parameters
        ----------
        arg1: string
            Suffix for output gtf file name
        arg2: tuple
            input Gtf files
        arg3: bool
            Overwrite if output file already exists.
        arg4: dict
            options to pass to stringtie
        """
        
        if len(args) < 1:
            print("ERROR: No input gtf for stringtie merge.")
            return ""
        
        #create path to output sam file
        fname=getFileBaseName(args[0])
        outDir=getFileDirectory(args[0])
        outGtfFile=os.path.join(outDir,fname+outFileSuffix+".gtf")
        
        if not overwrite:
            #check if file exists. return if yes
            if os.path.isfile(outGtfFile):
                print("The file "+outGtfFile+" already exists. Exiting..")
                return outGtfFile
        
        #Add merge flag, output file name and input bam
        newOpts={"--merge":"","-o":outGtfFile,"--":args}
        
        mergedOpts={**kwargs,**newOpts}
        
        #call stringtie
        status=self.runStringtie(verbose=verbose,quiet=quiet,logs=logs,objectid=objectid,**mergedOpts)
        
        if status:
            #check if sam file is present in the location directory of sraOb
            if checkFilesExists(outGtfFile):
                return outGtfFile
        else:
            return ""
        
        
        
            
    
    def runStringtie(self,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper for running stringtie
        
        Parameters
        ----------
        arg1: dict
            Options passed to stringtie
        
        Returns
        -------
        bool
            status of stringtie command.
        """
            
        #override existing arguments
        mergedArgsDict={**self.passedArgumentDict,**kwargs}
       
        stie_Cmd=['stringtie']
        #add options
        stie_Cmd.extend(parseUnixStyleArgs(self.validArgsList,mergedArgsDict))        
        
                
        #start ececution
        status=executeCommand(stie_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            printBoldRed("stringtie failed")
        
        #return status
        return status
    
    
    
class Cufflinks(Assembly):
    def __init__(self,referenceGTF="",**kwargs):
        """Stringtie constructor. Initialize stringtie parameters.
        
        Parameters
        ----------
        arg1: string
            Path to the reference gtf file.
        arg2: dict
            Options passed to stringtie. These could be overridden later when executing cufflinks.
        """
        super().__init__()
        self.programName="cufflinks"
        #check if stringtie exists
        if not checkDep([self.programName]):
            raise Exception("ERROR: "+ self.programName+" not found.")
            
        
        
        self.cufflinksArgsList=['-h','--help','-o','--output-dir','-p','--num-threads','--seed','-G','--GTF','-g','--GTF-guide','-M','--mask-file','-b','--frag-bias-correct','-u','--multi-read-correct','--library-type','--library-norm-method',
'-m','--frag-len-mean','-s','--frag-len-std-dev','--max-mle-iterations','--compatible-hits-norm','--total-hits-norm','--num-frag-count-draws','--num-frag-assign-draws','--max-frag-multihits','--no-effective-length-correction',
'--no-length-correction','-N','--upper-quartile-norm','--raw-mapped-norm','-L','--label','-F','--min-isoform-fraction','-j','--pre-mrna-fraction','-I','--max-intron-length','-a','--junc-alpha','-A','--small-anchor-fraction',
'--min-frags-per-transfrag','--overhang-tolerance','--max-bundle-length','--max-bundle-frags','--min-intron-length','--trim-3-avgcov-thresh','--trim-3-dropoff-frac','--max-multiread-fraction','--overlap-radius',
'--no-faux-reads','--3-overhang-tolerance','--intron-overhang-tolerance','-v','--verbose','-q','--quiet','--no-update-check']

        self.cuffcompareArgsList=['-h','-i','-r','-R','-Q','-M','-N','-s','-e','-d','-p','-C','-F','-G','-T','-V']
        self.cuffquantArgsList=['-o','--output-dir','-p','--num-threads','-M','--mask-file','-b','--frag-bias-correct','-u','--multi-read-correct','--library-type','-m','--frag-len-mean','-s','--frag-len-std-dev','-c','--min-alignment-count',
'--max-mle-iterations','-v','--verbose','-q','--quiet','--seed','--no-update-check','--max-bundle-frags','--max-frag-multihits','--no-effective-length-correction','--no-length-correction','--read-skip-fraction',
'--no-read-pairs','--trim-read-length','--no-scv-correction']
        self.cuffdiffArgsList=['-o','--output-dir','-L','--labels','--FDR','-M','--mask-file','-C','--contrast-file','-b','--frag-bias-correct','-u','--multi-read-correct','-p','--num-threads','--no-diff','--no-js-tests','-T','--time-series',
'--library-type','--dispersion-method','--library-norm-method','-m','--frag-len-mean','-s','--frag-len-std-dev','-c','--min-alignment-count','--max-mle-iterations','--compatible-hits-norm','--total-hits-norm',
' -v','--verbose','-q','--quiet','--seed','--no-update-check','--emit-count-tables','--max-bundle-frags','--num-frag-count-draws','--num-frag-assign-draws','--max-frag-multihits','--min-outlier-p','--min-reps-for-js-test',
'--no-effective-length-correction','--no-length-correction','-N','--upper-quartile-norm','--geometric-norm','--raw-mapped-norm','--poisson-dispersion','--read-skip-fraction','--no-read-pairs','--trim-read-length','--no-scv-correction']
        self.cuffnormArgsList=['-o','--output-dir','-L','--labels','--norm-standards-file','-p','--num-threads','--library-type','--library-norm-method','--output-format','--compatible-hits-norm','--total-hits-norm','-v','--verbose','-q','--quiet','--seed','--no-update-check']
        self.cuffmergeArgsList=['h','--help','-o','-g','–-ref-gtf','-p','–-num-threads','-s','-–ref-sequence']
        
        self.validArgsList=getListUnion(self.cufflinksArgsList,self.cuffcompareArgsList,self.cuffquantArgsList,self.cuffdiffArgsList,self.cuffnormArgsList,self.cuffmergeArgsList)
        
        #keep the passed arguments
        self.passedArgumentDict=kwargs
        
        #check the reference GTF
        if len(referenceGTF)>0 and checkFilesExists(referenceGTF):
            self.referenceGTF=referenceGTF
            self.passedArgumentDict['-g']=referenceGTF
    
    
    def performAssembly(self,inputBAM,outFileSuffix="_cufflinks",overwrite=True,**kwargs):
        """Function to run cufflinks with BAM file as input.
                
        Parameters
        ----------
        arg1: string
            path to bam file
        arg2: string
            Suffix for the output gtf file
        arg3: bool
            Overwrite if output file already exists.
        arg4: dict
            Options to pass to stringtie. This will override the existing options self.passedArgumentDict (only replace existing arguments and not replace all the arguments).
            
        Returns
        -------
        string
            path to output GTF file
        
        """
        
        #create path to output file
        fname=getFileBaseName(inputBAM)
        outDir=getFileDirectory(inputBAM)
        outGtfFile=os.path.join(outDir,fname+outFileSuffix+".gtf")
        
        """
        Handle overwrite
        """
        if not overwrite:
            #check if file exists. return if yes
            if os.path.isfile(outGtfFile):
                print("The file "+outGtfFile+" already exists. Exiting..")
                return outGtfFile
            
        #Add output file name and input bam
        newOpts={"-o":outDir,"--":(inputBAM,)}
        mergedOpts={**kwargs,**newOpts}
        
        #call cufflinks
        status=self.runCufflinks(**mergedOpts)
        
        if status:
            #move outDir/transcripts.gtf to outfile
            moveFile(os.path.join(outDir,"transcripts.gtf"),outGtfFile)
            #check if sam file is present in the location directory of sraOb
            if checkFilesExists(outGtfFile):
                return outGtfFile
        else:
            return ""
    
    def runCuffCommand(self,command,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper for running cuff* commands
        
        Parameters
        ----------
        command: string
            the command name
        arg2: dict
            Options passed to cuff command
        
        Returns
        -------
        bool
            return status of the command.
        """
        validCommands=['cuffcompare','cuffdiff', 'cufflinks', 'cuffmerge', 'cuffnorm', 'cuffquant']
        if command in validCommands:
            #override existing arguments
            mergedArgsDict={**self.passedArgumentDict,**kwargs}
       
            cuff_Cmd=[command]
            #add options
            cuff_Cmd.extend(parseUnixStyleArgs(self.validArgsList,mergedArgsDict))        
                  
            #start ececution
            status=executeCommand(cuff_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
            if not status:
                printBoldRed("cufflinks failed")
                #return status
            return status
        else:
            printBoldRed("Unknown command {}"+command)
            return False
    
    
    def runCufflinks(self,verbose=False,quiet=False,logs=True,objectid="NA",**kwargs):
        """Wrapper for running cufflinks
        
        Parameters
        ----------
        arg1: dict
            Options passed to cufflinks
        
        Returns
        -------
        bool
            status of cufflinks command.
        """
            
        #override existing arguments
        mergedArgsDict={**self.passedArgumentDict,**kwargs}
       
        cufflinks_Cmd=['cufflinks']
        #add options
        cufflinks_Cmd.extend(parseUnixStyleArgs(self.validArgsList,mergedArgsDict))        
        
        
        #start ececution
        status=executeCommand(cufflinks_Cmd,verbose=verbose,quiet=quiet,logs=logs,objectid=objectid)
        if not status:
            printBoldRed("cufflinks failed")
        #return status
        return status