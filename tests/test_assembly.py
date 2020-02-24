#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 16:24:44 2019

@author: usingh
"""
import pytest
from pyrpipe import assembly
from pyrpipe import pyrpipe_utils as pu
from testingEnvironment import testSpecs
import os

testVars=testSpecs()

def test_assembly():
    #create assemble object
    aob=assembly.Assembly()
    assert aob.category == "Assembler","Failed Assembly init"

def test_stringtie():
    bam=testVars.hisatSortedBam
    gtf=testVars.gtf
    stie=assembly.Stringtie(reference_gtf=gtf)
    result=stie.perform_assembly(bam,out_dir=testVars.testDir, objectid="test")
    assert pu.check_files_exist(result)==True, "Failed stringtie"
    
def test_cufflinks():
    bam=testVars.hisatSortedBam
    gtf=testVars.gtf
    cl=assembly.Cufflinks(reference_gtf=gtf)
    result=cl.perform_assembly(bam,out_dir=testVars.testDir, objectid="test")
    assert pu.check_files_exist(result)==True, "Failed stringtie"
    
    
def test_trinityBam():
    tr=assembly.Trinity()
    #bam=testVars.hisatSortedBam
    bam=testVars.portcullisBam
    tr_opts={"--genome_guided_bam":bam,
             "--genome_guided_max_intron":"10000",
            "--output":testVars.testDir+"/trinity_testoutbam",
            "--max_memory": "2G"
            }
    st=tr.run_trinity(verbose=True,**tr_opts)
    assert st==True, "Failed trinity with bam"

def test_trinityFQ():
    tr=assembly.Trinity()
    tr_opts={"--seqType":"fq","--left":testVars.fq1,
            "--right":testVars.fq2,
            "--output":testVars.testDir+"/trinity_testoutfq",
            "--max_memory":"3G"}
    st=tr.run_trinity(**tr_opts)
    assert st==True, "Failed trinity with fq"


