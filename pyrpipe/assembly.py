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
        arg3: dict
            Options to pass to stringtie. This will override the existing options self.passedArgumentDict (only replace existing arguments and not replace all the arguments).
            
        Returns
        -------
        string
            path to output GTF file
        """
        
        #create path to output sam file
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
        newOpts={"-o":outGtfFile,"--":inputBAM}
        mergedOpts={**kwargs,**newOpts}
        print("MS:"+str(mergedOpts))
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
        log=""
        try:
            for output in executeCommand(stie_Cmd):
                #print (output)    
                log=log+str(output)
            #save to a log file
            
        except subprocess.CalledProcessError as e:
            print ("Error in command...\n"+str(e))
            #save error to error.log file
            return False        
        #return status
        return True
        
        
    #function to run stringtie
    def runStringtieold(self,inputBAM,proc,referenceGTF=False,outDir=False,deleteInputBam=False,outSuffix="_stie"):
                
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