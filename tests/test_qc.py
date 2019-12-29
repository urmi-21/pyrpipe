#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 28 01:18:25 2019

@author: usingh
"""

from pyrpipe import qc
from pyrpipe import pyrpipe_utils as pu
from testingEnvironment import testSpecs
import os

testVars=testSpecs()
fq1=testVars.fq1
fq2=testVars.fq2
rRNAfasta=testVars.rRNAfa
    
def test_RNASeqQC():
    ob=qc.RNASeqQC()
    assert ob.category=="RNASeqQC", "RNASeqQC failed"
    
    
def test_trimgalore():
    tg=qc.Trimgalore()
    
    #run tg
    tgOpts={"--cores": "10", "-o":testVars.testDir, "--paired":"", "--": (fq1,fq2)}
	     #remove fastq files}
    st=tg.run_trimgalore(**tgOpts)
    assert st==True, "Trimgalore failed"
    
def test_bbduk():
    bd=qc.BBmap()
    bdOpts={"in":fq1,"in2":fq2,"out":testVars.testDir+"/bdo1.fq","out2":testVars.testDir+"/bdo2.fq"}
    st=bd.run_bbduk(**bdOpts)
    assert st==True, "BBDUK failed"
    

def test_bbsplit():
    bs=qc.BBmap()
    #build index
    bsOpts={"ref_x":rRNAfasta,"path": testVars.testDir}
    st=bs.run_bbsplit(**bsOpts)
    assert st==True, "BBSplit index failes"
    
    #run split
    bsOpts={"in1":fq1,"in2":fq2,"outu1":testVars.testDir+"/bso1.fq","outu2":testVars.testDir+"/bso2.fq","path":testVars.testDir}
    st==bs.run_bbsplit(**bsOpts)
    assert st==True, "BBSplit failed"