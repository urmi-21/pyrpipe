#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 16:24:44 2019

@author: usingh
"""

from pyrpipe import sra
from testingEnvironment import testSpecs
import os

testVars=testSpecs()

def test_sra():
    srrID=testVars.srr
    print ("Testing SRA...")
    testDir=testVars.testDir
    newOb=sra.SRA(srrID,testDir)
    assert newOb.srr_accession==srrID,"Failed SRA init"
    assert newOb.download_sra()==True, "Failed to download SRA file"
    assert newOb.sraFileExistsLocally()==True, "Failed to locate .sra file on disk"
    ##test fasterq-dump
    assert newOb.run_fasterqdump(delete_sra=True,**{"-f":"","-t":testDir})==True, "Failed FQ dump"
    assert newOb.fastqFilesExistsLocally()==True, "Failed to locate .fastq files on disk"
    assert newOb.sraFileExistsLocally()!=True, "Failed to delete .sra files from disk"
    #delete downloaded files
    assert newOb.delete_fastq()==True, "Failed to delete .fastq files from disk"