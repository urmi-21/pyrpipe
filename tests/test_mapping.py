#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 23:32:09 2019

@author: usingh
"""

import pytest
from pyrpipe import mapping
from pyrpipe import pyrpipe_utils as pu
from testingEnvironment import testSpecs
import os

testVars=testSpecs()

def test_mapping():
    #create assemble object
    mob=mapping.Aligner()
    assert mob.category == "Aligner","Failed Mapping init"
    
def test_hisat2():
    #test hisat build and hisat mapping
    #fq1="tests/test_files/athaliana/fastq1/hisat2_sorted.bam"
    #gtf="tests/test_files/athaliana/genome/Arabidopsis_thaliana.TAIR10.45_1and2.gtf"
    #stie=assembly.Stringtie(reference_gtf=gtf)
    #result=stie.perform_assembly(bam,out_dir=testVars.testDir, objectid="test")
    #assert pu.check_files_exist(result)==True, "Failed stringtie"
    assert False==True, "Failed stringtie"