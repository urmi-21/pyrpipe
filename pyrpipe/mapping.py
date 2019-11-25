#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 19:53:42 2019

@author: usingh
contains classes of RNA-Seq mapping programs
"""

from myutils import *

class HISAT2:
    def __init__(self,hisat2Index="",*args):
        """HISAT2 constructor. Initialize hisat2's index and other parameters.
        Parameters
        ----------
        
        """
        
        
        #check index exists
        if len(hisat2Index)>0 and checkHisatIndex(hisat2Index):
            print("Found HISAT2 index files.")
            self.hisat2Index=hisat2Index
        elif args :
            print("Building HISAT2 index with "+ str(args))
            indexPath,indexName=os.path.split(hisat2Index)
            self.buildHisat2Index(indexPath,indexName,*args)
        else:
            raise Exception("ERROR: Invalid HISAT2 index. Please run build index to generate an index.")
            
    def buildHisat2Index(self,indexPath,indexName,*args):
        print("Building hisat index...")
        
        
    
    def runHisat2(self,sraOb,**kwargs):
        """Run HISAT2 using and SRA object and produce .bam file as result. The HISAT2 index used will be self.hisat2Index.
        All output will be written to SRA.location by default.
        
        Parameters
        ----------
        arg1: SRA object
            An object of type SRA which contains valid fastq files
        
        arg2: dict
            arguments to be passed to hisat2
        """
            





class STAR:
    def __init__(self,starIndex):
        """STAR constructor. Initialize star's index and other parameters.
        """

if __name__ == "__main__":
    #test
    
    hs=HISAT2("ssa","as","dsa","","dsadrr")
    
    print ("done")