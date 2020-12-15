#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 11:18:11 2020

@author: usingh
"""

from pyrpipe import sra,pyrpipe_session
from testingEnvironment import testSpecs

testVars=testSpecs()

def test_sra():
    srrID=testVars.srr
    print ("Testing SRA...")
    testDir=testVars.testDir
    newOb=sra.SRA(srrID,testDir)
    assert newOb.srr_accession==srrID,"Failed SRA init"
    assert newOb.fastq_exists()==True, "Failed to locate .fastq files on disk"
    #delete fastq
    assert newOb.delete_fastq()==True, "failed to delete fastq"
    
    
    assert newOb.download_sra()==True, "Failed to download SRA file"
    assert newOb.sra_exists()==True, "Failed to locate .sra file on disk"
    ##test fasterq-dump
    assert newOb.run_fasterqdump(delete_sra=True,**{"-f":"","-t":testDir})==True, "Failed FQ dump"
    assert newOb.sra_exists()==False, "Failed to delete .sra files from disk"
    
    st=pyrpipe_session.save_session(filename="mySession",add_timestamp=False,out_dir=testVars.testDir)
    assert st==True, "Session save failed"
    st=pyrpipe_session.restore_session(testVars.testDir+"/mySession.pyrpipe")
    assert st==True, "Session restore failed"
    
    assert newOb.fastq_exists()==True, "Failed to locate .fastq files on disk after restore"
    assert newOb.delete_fastq()==True, "failed to delete fastq"
    
