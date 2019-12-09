#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:21:01 2019

@author: usingh
"""

from pyrpipe.pyrpipe_utils import *

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
        """Function to run stringtie with an object of SRA class.
                
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
        
    def performStringtieMerge(self,*args,outFileSuffix="_stringtieMerge",overwrite=True,**kwargs):
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
        status=self.runStringtie(**mergedOpts)
        
        if status:
            #check if sam file is present in the location directory of sraOb
            if checkFilesExists(outGtfFile):
                return outGtfFile
        else:
            return ""
        
        
        
            
    
    def runStringtie(self,**kwargs):
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
        
        print("Executing:"+" ".join(stie_Cmd))
        
        #start ececution
        status=executeCommand(stie_Cmd)
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
            Options passed to stringtie. These could be overridden later when executing stringtie.
        """
        super().__init__()
        self.programName="cufflinks"
        #check if stringtie exists
        if not checkDep([self.programName]):
            raise Exception("ERROR: "+ self.programName+" not found.")
        self.validArgsList=[]
        
        #keep the passed arguments
        self.passedArgumentDict=kwargs
        
        #check the reference GTF
        if len(referenceGTF)>0 and checkFilesExists(referenceGTF):
            self.referenceGTF=referenceGTF
            self.passedArgumentDict['-g']=referenceGTF
    
    
    def performAssembly(self,inputBAM,outFileSuffix="_cufflinks",overwrite=True,**kwargs):
        """
        """
        pass
    
    
    
    def runCufflinks(self,**kwargs):
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
        
        print("Executing:"+" ".join(cufflinks_Cmd))
        
        #start ececution
        status=executeCommand(cufflinks_Cmd)
        if not status:
            printBoldRed("cufflinks failed")
        #return status
        return status