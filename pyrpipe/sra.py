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
        The **kwargs is for passing arguments to the prefetch command
        Parameters
        ----------
        arg1 : string
            location on disk to downloaded the file. prefetch will create a folder with name same as <srrAccession> under the path specified.
            Default is the current directory
        
        Returns
        -------
        bool
            Return status of the prefetch command. 0 for successfull download and 1 for failiure

        Examples
        --------
        >>> object.downloadSRAFile()
        True
        """
        
        print("Downloading "+self.srrAccession+" ...")
        
        #scan for prefetch arguments
        prefetchArgsList=['-f','-t','-l','-n','-s','-R','-N','-X','-o','-a','--ascp-options','-p','--eliminate-quals','-c','-o','-O','-h','-V','-L','-v','-q']
        pathFound=False
        
        #append srrid at the end of path
        if '-O' in kwargs:
            print("Ignoring -O flag."+" location is: "+self.location)
            #delete -O parameter
            del kwargs['-O']
            '''
            pathFound=True
            value=kwargs['-O']
            value=os.path.join(value,self.srrAccession)
            kwargs['-O']=value
            sraLocation=value
            '''

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
        return True
           
    
    def sraFileExistsLocally(self):
        if hasattr(self,'localSRAPath'):
            return(os.path.isfile(self.localSRAFilePath))
        else:
            return False
    
        

if __name__ == "__main__":
    #test
    newOb=SRA('SRR10408795',"/home/usingh/work/urmi/hoap/test")
    print(newOb.getSrrAccession())
    
    print(newOb.sraFileExistsLocally())
    newOb.downloadSRAFile(**{"-O": "/home/usingh/work/urmi/hoap", "Attr2": "Val2","-q":""})
    print(newOb.sraFileExistsLocally())
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    