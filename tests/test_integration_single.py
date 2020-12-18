#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 17:28:58 2020

@author: usingh
"""

from pyrpipe import sra,qc,mapping,assembly
from testingEnvironment import testSpecs
import os

testVars=testSpecs()
fq1=testVars.fq1
fq2=testVars.fq2
rRNAfasta=testVars.rRNAfa

srr='ERR3770564' #single end arabidopsis data

workingDir=testVars.testDir


#create objects
bbdOpts={"ktrim":"r","k":"23","mink":"11","qtrim":"'rl'","trimq":"10","ref":testVars.bbdukAdapters}
bbdOb=qc.BBmap(None,**bbdOpts)
tg=qc.Trimgalore()
bt=mapping.Bowtie2(index=testVars.testDir+"/btIndex",genome=testVars.genome)
hsOpts={"--dta-cufflinks":"","-p":"8"}
hs=mapping.Hisat2(index=testVars.testDir+"/hisatindex",genome=testVars.genome,**hsOpts)
star=mapping.Star(index=os.path.join(testVars.testDir,"starIndex"),genome=testVars.genome)
stie=assembly.Stringtie()


#sra ob 
sraOb=sra.SRA(srr,workingDir)
st=sraOb.fastq_exists()
assert st==True,"fasterq-dump failed"


def test_pipeline1():    
    st=sraOb.trim(bbdOb).align(hs).assemble(stie)
    assert st!=None,"pipeline 1 failed"
    
def test_pipeline2():    
    st=sraOb.trim(tg).align(hs).assemble(stie)  
    assert st!=None,"pipeline 1 failed"

  
def test_pipeline3():    
    st=sraOb.trim(tg).align(star).assemble(stie)
    assert st!=None,"pipeline 1 failed"

    
def test_pipeline4():    
    st=sraOb.align(star)
    assert st!=None,"pipeline 1 failed"
