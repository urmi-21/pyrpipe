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

#mergedBam=""


def test_samtools():
    #test sam to sorted bam
    sm=tools.Samtools()    
    sortedBam=sm.sam_sorted_bam(testVars.hisatSam,out_dir=testVars.testDir)
    print("check:"+sortedBam)
    st=pu.check_files_exist(sortedBam)
    assert st==True, "Failed to convert sam to sorted bam"
    
    #test merge
    mergedBam=sm.merge_bam(testVars.hisatSortedBam,testVars.starSortedBam,out_dir=testVars.testDir,**{"-f":""})
    st=pu.check_files_exist(mergedBam)
    assert st==True, "Failed to merge bam"
    
def test_portcullis():
    pc=tools.Portcullis()
    port_out=pc.run_portcullisFull(testVars.genome,testVars.hisatSortedBam,out_dir=testVars.testDir)
    st=pu.check_paths_exist(port_out)
    assert st==True, "Failed portcullis run"
    
    
def test_mikado():
    gtfdir=testVars.mikadofiles
    out_dir=testVars.testDir+"/mikadoout"
    mk=tools.Mikado()
    
    listfile=mk.createMikadoGTFlist("mikadolist",out_dir,gtfdir)
    st=pu.check_files_exist(listfile)
    assert st==True, "Mikado list failed"