#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 16:24:44 2019

@author: usingh
"""

from pyrpipe import assembly
from pyrpipe import pyrpipe_utils as pu
from testingEnvironment import testSpecs

testVars=testSpecs()

def test_assembly():
    #create assemble object
    aob=assembly.Assembly()
    assert aob._category == "Assembler","Failed Assembly init"

def test_stringtie():
    bam=testVars.hisatSortedBam
    gtf=testVars.gtf
    stie=assembly.Stringtie()
    result=stie.perform_assembly(bam,out_dir=testVars.testDir, objectid="test",reference_gtf=gtf)
    assert pu.check_files_exist(result)==True, "Failed stringtie"
    
def test_cufflinks():
    bam=testVars.hisatSortedBam
    gtf=testVars.gtf
    cl=assembly.Cufflinks()
    result=cl.perform_assembly(bam,out_dir=testVars.testDir, objectid="test",reference_gtf=gtf)
    assert pu.check_files_exist(result)==True, "Failed cufflinks"
    
    



