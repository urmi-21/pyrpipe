#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 15:56:20 2019

@author: usingh
"""

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