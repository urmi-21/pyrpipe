#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 15:45:26 2019

@author: usingh


"""

from myutils import *
import os

class SRA:
    def __init__(self,srrAccession):
        self.srrAccession=srrAccession
        
    def getSrrAccession(self):
            return self.srrAccession
    
  
    def downloadSRAFile(self,**kwargs):
        """This function downloads .sra file from NCBI SRA servers using the prefetch command.

        NCBI sra-toolkit 2.9 or higher must be installed on the system in order to use prefetch. 
        prefetch will create a folder with name same as <srrAccession> under the path specified.
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
        sraLocation=""
        
        #append srrid at the end of path
        if '-O' in kwargs:
            pathFound=True
            value=kwargs['-O']
            value=os.path.join(value,self.srrAccession)
            kwargs['-O']=value
            sraLocation=value
            
        
        
        prefetch_Cmd=['prefetch']
        prefetch_Cmd.extend(parseUnixStyleArgs(prefetchArgsList,kwargs))
        #if path is not specified, used current directory
        if not pathFound:
            prefetch_Cmd.extend(['-O',os.path.join(os.getcwd(),self.srrAccession)])
            sraLocation=os.path.join(os.getcwd(),self.srrAccession)
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
        self.localSRAPath=os.path.join(sraLocation,self.srrAccession+".sra")
        #validate path exists
        if not checkFilesExists(self.localSRAPath):
            printBoldRed("Error downloading file. File "+self.localSRAPath+" does not exist!!!")
            return False
        
        print ("Downloaded file: "+self.localSRAPath+" {0} ".format(getFileSize(self.localSRAPath)))
        return True
    
    def sraFileExistsLocally(self):
        if hasattr(self,'localSRAPath'):
            return(os.path.isfile(self.localSRAPath))
        else:
            return False
        
        

if __name__ == "__main__":
    #test
    newOb=SRA('SRR10408795')
    print(newOb.getSrrAccession())
    
    print(newOb.sraFileExistsLocally())
    newOb.downloadSRAFile(**{"-O": "/home/usingh/work/urmi/hoap/test", "Attr2": "Val2","-q":""})
    print(newOb.sraFileExistsLocally())
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    