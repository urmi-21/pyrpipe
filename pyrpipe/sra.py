#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 15:45:26 2019

@author: usingh


"""

class SRA:
    def __init__(self,srrAccession):
        self.srrAccession=srrAccession
        
    def getSrrAccession(self):
            return self.srrAccession
    
  
    def downloadSRAFile(self,**kwargs):
        """This function downloads .sra file from NCBI SRA servers using the prefetch command.

        NCBI sra-toolkit must be installed on the system in order to use prefetch. 
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
        
        
        

if __name__ == "__main__":
    #test
    newOb=SRA('SRR002328')
    print(newOb.getSrrAccession())
    newOb.downloadSRAFile(name="yasoobd")