#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 15:45:26 2019

@author: usingh


"""

from myutils import *
import os

class SRA:
    def __init__(self,srrAccession,location=os.getcwd()):
        """Init SRA object.
        Parameters
        ----------
        arg1 : string
            A valid SRR accession id
        arg2 : String
            Path where all data related to this object (e.g. .sra files, metadata, fastq files) will be stored. 
            Default value of the path will be "./<SRR_accession>". <SRR_accession> is added at the end of the path
            so that final location is location/<SRR_accession>.
            For consistency, location and SRR Accession id are not allowed to be modified.
        
        """
        self.srrAccession=srrAccession
        #append the SRR accession to the location
        self.location=os.path.join(location,self.srrAccession)
        
    def __setattr__(self, name, value):
        if name not in self.__dict__:
            self.__dict__[name] = value
        else:
            raise Exception("Can not modify this field")
        
    def getSrrAccession(self):
            return self.srrAccession
    
  
    def downloadSRAFile(self,**kwargs):
        """This function downloads .sra file from NCBI SRA servers using the prefetch command.

        NCBI sra-toolkit 2.9 or higher must be installed on the system in order to use prefetch. 
        prefetch will create a folder with name same as <srrAccession> under the location (path) specified.
        The path of downloaded file is saved in the object as localSRAPath. This localSRAPath is then used
        by other functions to access the downloaded data. 
        The **kwargs is for passing arguments to the prefetch command.
        
        Parameters
        ----------
        arg1:
            dict containing prefetch arguments
        
        Returns
        -------
        bool
            Return status of the prefetch command. True if successful download and False if failed.

        Examples
        --------
        >>> object.downloadSRAFile()
        True
        """
        
        print("Downloading "+self.srrAccession+" ...")
        
        #scan for prefetch arguments
        prefetchArgsList=['-f','-t','-l','-n','-s','-R','-N','-X','-o','-a','--ascp-options','-p','--eliminate-quals','-c','-o','-O','-h','-V','-L','-v','-q']
        pathFound=False
        
        #ignore location and file name arguments if given
        if '-O' in kwargs:
            print("Ignoring -O flag."+" location is: "+self.location)
            #delete -O parameter
            del kwargs['-O']
        if '-o' in kwargs:
            print("Ignoring -o flag."+" File name is: "+self.srrAccession)
            #delete -o parameter
            del kwargs['-o']
            

        prefetch_Cmd=['prefetch']
        prefetch_Cmd.extend(parseUnixStyleArgs(prefetchArgsList,kwargs))
        prefetch_Cmd.extend(['-O',self.location])
        prefetch_Cmd.append(self.srrAccession)
        print("Executing:"+" ".join(prefetch_Cmd))
        
        
        log=""
        try:
            for output in executeCommand(prefetch_Cmd):
                print (output)    
                log=log+str(output)
            #save to a log file

        except subprocess.CalledProcessError as e:
            print ("Error in command...\n"+str(e))
            #save error to error.log file
            return False
        
        #store path to the downloaded sra file
        self.localSRAFilePath=os.path.join(self.location,self.srrAccession+".sra")
        #validate path exists
        if not checkFilesExists(self.localSRAFilePath):
            printBoldRed("Error downloading file. File "+self.localSRAFilePath+" does not exist!!!")
            return False
        
        print ("Downloaded file: "+self.localSRAFilePath+" {0} ".format(getFileSize(self.localSRAFilePath)))
        #test if file is paired or single end
        if isPairedSRA(self.localSRAFilePath):
            self.layout="PAIRED"
        else:
            self.layout="SINGLE"
            
        
        return True
    
    
           
    
    def sraFileExistsLocally(self):
        if hasattr(self,'localSRAFilePath'):
            return(os.path.isfile(self.localSRAFilePath))
        else:
            return False
    
    def runFasterQDump(self,deleteSRA=False,**kwargs):
        """Execute fasterq-dump to convert .sra file to fastq files.
        The fastq files will be stored in the same directory as the sra file.
        
        Parameters
        ----------
        arg1: bool
            delete sra file after completion
            
        arg2: dict
            A dict containing fasterq-dump arguments
        
        Returns
        -------
        bool
            Return status of the fasterq-dump command. True if successful download and False if failed.

        Examples
        --------
        >>> object.runFasterQDump()
        True
        """
        
        #first check is sra exists
        if not self.sraFileExistsLocally():
            print ("Error executing fasterq-dump: .sra file not found. Please run downloadSRAFile().")
            return False
        #else directly run fasterq-dump on accession ?
        
        fasterqdumpArgsList=['-f','-t','-s','-N','-X','-a','-p','-c','-o','-O','-h','-V','-L','-v','-q','-b','-m','-e','-x','-S','-3','-P','-M','-B','--option-file','--strict','--table','--include-technical','--skip-technical','--concatenate-reads']
        
        
        
        #ignore location and file name arguments if given
        if '-O' in kwargs:
            print("Ignoring -O flag."+" location is: "+self.location)
            #delete -O parameter
            del kwargs['-O']
        if '-o' in kwargs:
            print("Ignoring -o flag."+" File name is: "+self.srrAccession)
            #delete -o parameter
            del kwargs['-o']
        
        
        #execute command
        
        fstrqd_Cmd=['fasterq-dump']
        fstrqd_Cmd.extend(parseUnixStyleArgs(fasterqdumpArgsList,kwargs))
        #add location
        fstrqd_Cmd.extend(['-O',self.location])
        #add output filename. output will be <srrAccession>.fastq or <srrAccession>_1.fastq and <srrAccession>_2.fastq
        fstrqd_Cmd.extend(['-o',self.srrAccession+".fastq"])
        fstrqd_Cmd.append(self.localSRAFilePath)
        print("Executing:"+" ".join(fstrqd_Cmd))
        
        log=""
        try:
            for output in executeCommand(fstrqd_Cmd):
                print (output)    
                log=log+str(output)
            #save to a log file

        except subprocess.CalledProcessError as e:
            print ("Error in command...\n"+str(e))
            #save error to error.log file
            return False
        
        #check if fastq files are downloaded 
        if(self.layout=="SINGLE"):
            self.localfastqPath=os.path.join(self.location,self.srrAccession+".fastq")
            
            if not checkFilesExists(self.localfastqPath):
                printBoldRed("Error running fasterq-dump file. File "+self.localfastqPath+" does not exist!!!")
                return False
            
        else:
            self.localfastq1Path=os.path.join(self.location,self.srrAccession+"_1.fastq")
            self.localfastq2Path=os.path.join(self.location,self.srrAccession+"_2.fastq")
            
            if not checkFilesExists(self.localfastq1Path,self.localfastq2Path):
                printBoldRed("Error running fasterq-dump file. File "+self.localfastq1Path+" does not exist!!!")
                return False
        
        return True
        
        

if __name__ == "__main__":
    #test
    testDir="/home/usingh/work/urmi/hoap/test"
    newOb=SRA('SRR10408795',testDir)
    print(newOb.sraFileExistsLocally())
    newOb.downloadSRAFile(**{})
    print(newOb.sraFileExistsLocally())
    #-e numthreads, -S split files, -t temp directory
    newOb.runFasterQDump(**{"-e":"8","-S":"","--skip-technical":"","-t":testDir})
    
    
    sraOb2=SRA("ERR3520221",testDir)
    sraOb2.downloadSRAFile()
    #-e numthreads, -S split files, -t temp directory
    sraOb2.runFasterQDump(**{"-e":"8","-S":"","--skip-technical":"","-t":testDir})
    
    import mapping
    hs=mapping.HISAT2("xyz","a")
    
    print ("done")
    
    
    
    
    
    
    
    
    
    
    