#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:21:01 2019

@author: usingh
"""

from pyrpipe.myutils import *

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
            Path to the reference gtf file. Note: -G is not in the validArgsList
        arg2: dict
            Options passed to stringtie. These could be overridden later when executing stringtie.
        """
        super().__init__()
        self.programName="stringtie"
        #check if stringtie exists
        if not checkDep([self.programName]):
            raise Exception("ERROR: "+ self.programName+" not found.")
        self.validArgsList=['--version','--conservative','--rf','--fr','-o','-l',
                            '-f','-L','-m','-a','-j','-t','-c','-s','-v','-g','-M',
                            '-p','-A','-B','-b','-e','-x','-u','-h','--merge','-F','-T','-i']
        
        #keep the passed arguments
        self.passedArgumentDict=kwargs
        
    def performAssembly():
        """Function to run stringtie with an object of SRA class.
        """
        pass
            
 #function to run stringtie
    def runStringtie(self,inputBAM,proc,referenceGTF=False,outDir=False,deleteInputBam=False,outSuffix="_stie"):
                
        """A wrapper to run stringtie. All arguments are passed through **kwargs.
        This function is called by performAssembly() or could be called independently.
        
        Returns
        -------
        string
            Returns the absolute path og the GTF file if the assembly was successful. Returns empty string otherwise.
        """
        
        #check if bam file exists
        if not checkFilesExists(inputBAM):
            raise Exception("ERROR: Can not find bam file: "+ inputBAM)
        
        #use bam file directory as default directory
        if not outDir:
            outDir=os.path.split(inputBAM)[0]
        fname=os.path.split(inputBAM)[1]
        outGtfFile=os.path.join(outDir,fname+outSuffix+'.gtf')
        
        #stie command
        stie_Cmd=['stringtie',inputBAM,'-p',str(proc)]
        
        if referenceGTF and checkFilesExists(referenceGTF):
            stie_Cmd.extend(['-G',referenceGTF])
                
        stie_Cmd.extend(['-o',outGtfFile])
        print("Executing: "+" ".join(stie_Cmd))
        
        try:
            for output in executeCommand(stie_Cmd):
                print (output)
        except subprocess.CalledProcessError as e:
            print ("Error in command")
            return ""
        
        
        if deleteInputBam:
            if not deleteFileFromDisk(inputBAM):
                return ""
        
        #check gtf file
        if not checkFilesExists(outGtfFile):
            return ""
        
        
        
        return outGtfFile