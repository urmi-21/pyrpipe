#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 23:32:09 2019

@author: usingh
"""

from pyrpipe import mapping
from testingEnvironment import testSpecs
import os

testVars=testSpecs()

def test_mapping():
    #create assemble object
    mob=mapping.Aligner()
    assert mob._category == "Aligner","Failed Mapping init"
    

def test_hisat2():
    #test hisat build and hisat mapping
    hs=mapping.Hisat2(index=os.path.join(testVars.testDir,"hisatindex"),genome=testVars.genome)
    assert hs.check_index()==True, "Failed hisat buildindex"
    #perform alignment without sraobject
    kwargs={"-1":testVars.fq1,"-2":testVars.fq2,"-S":testVars.testDir+"/hisatTest.sam","--dta-cufflinks":"","-p":"4"}
    st=hs.run(**kwargs)
    assert st==True, "Failed to run hisat"


def test_star():
    star=mapping.Star(index=os.path.join(testVars.testDir,"starIndex"),genome=testVars.genome)
    assert star.check_index()==True, "Failed star build_index"
    
    #perform alignment without sraobject
    #create outdir
    outdir=os.path.join(testVars.testDir,"starout")
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    opts={"--outFilterType":"BySJout",
            "--runThreadN":"8",
            "--outSAMtype": "BAM SortedByCoordinate",
            "--outFileNamePrefix":outdir+"/",
            "--readFilesIn":testVars.fq1+" "+testVars.fq2
            }
    st=star.run(**opts)
    assert st==True, "Failed to run star"



def test_bowtie():
    bt=mapping.Bowtie2(index=os.path.join(testVars.testDir,"btIndex"),genome=testVars.genome)
    assert bt.check_index()==True, "Failed bowtie2 build_index"
    
    opts={"-1":testVars.fq1,"-2":testVars.fq2,"-S":testVars.testDir+"/bowtie2.sam"}
    st=bt.run(**opts)
    assert st==True, "Failed to run bowtie2"

