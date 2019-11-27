#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:21:01 2019

@author: usingh
"""

from pyrpipe.myutils import *

class Stringtie:
    def __init__(self):
        self.programName="stringtie"
        #check if hisat2 exists
        if not checkDep([self.programName]):
            raise Exception("ERROR: "+ self.programName+" not found.")
            
            
 #function to run stringtie
    def runStringtie(self,inputBAM,proc,referenceGTF=False,outDir=False,deleteInputBam=False,outSuffix="_stie"):
        
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
            return False,"NA"
        
        
        if deleteInputBam:
            if not deleteFileFromDisk(inputBAM):
                return False,""
        
        #check gtf file
        if not checkFilesExists(outGtfFile):
            return False,"NA"
        
        
        
        return True,outGtfFile