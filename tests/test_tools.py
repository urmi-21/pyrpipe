#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 16:40:00 2019

@author: usingh
"""

from pyrpipe import tools
from pyrpipe import pyrpipe_utils as pu
from testingEnvironment import testSpecs
#import os

testVars=testSpecs()


def test_samtools():
    #test sam to sorted bam
    sm=tools.Samtools(**{"-@":"2"})    
    sortedBam=sm.sam_sorted_bam(testVars.hisatSam,out_dir=testVars.testDir)
    print("check:"+sortedBam)
    st=pu.check_files_exist(sortedBam)
    assert st==True, "Failed to convert sam to sorted bam"
    
    #test merge
    bam_list=[testVars.hisatSortedBam,testVars.starSortedBam]
    mergedBam=sm.merge_bam(bam_list,out_dir=testVars.testDir)
    st=pu.check_files_exist(mergedBam)
    assert st==True, "Failed to merge bam"
 


    
    
    
    
    
