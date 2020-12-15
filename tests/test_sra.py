#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 16:24:44 2019

@author: usingh
"""

from pyrpipe import sra
from testingEnvironment import testSpecs

testVars=testSpecs()

def test_sra():
    srrID=testVars.srr
    print ("Testing SRA...")
    testDir=testVars.testDir
    newOb=sra.SRA(srrID,testDir+"/testsra")
    assert newOb.srr_accession==srrID,"Failed SRA init"
    assert newOb.download_sra()==True, "Failed to download SRA file"
    assert newOb.sra_exists()==True, "Failed to locate .sra file on disk"
    ##test fasterq-dump
    assert newOb.run_fasterqdump(delete_sra=True,**{"-f":"","-t":testDir})==True, "Failed FQ dump"
    assert newOb.fastq_exists()==True, "Failed to locate .fastq files on disk"
    assert newOb.delete_sra()==True, "Failed to delete .sra files from disk"    
    assert newOb.sra_exists()==False, "Deleted SRA file exists. "
    #delete downloaded files
    assert newOb.delete_fastq()==True, "Failed to delete .fastq files from disk"