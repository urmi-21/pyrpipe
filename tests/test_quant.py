#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 18:26:06 2020

@author: usingh
"""
from pyrpipe import quant
from testingEnvironment import testSpecs
#import os

testVars=testSpecs()

def test_mapping():
    #create assemble object
    qob=quant.Quant()
    assert qob._category == "Quantification","Failed Quant init"
 
def test_kallisto():
    kl=quant.Kallisto(index=testVars.testDir+"/kallistoIndex/kalIndex",transcriptome=testVars.cdna)
    assert kl.check_index()==True, "Failed kallisto build_index"
    args=(testVars.fq1,testVars.fq2)
    opts={"-o":testVars.testDir+"/kalOut"}
    st=kl.run(*args,subcommand='quant',target='tests/testout/kalOut/abundance.tsv',**opts)
    assert st==True, "Failed to run kallisto"
    
    
def test_salmon():
    sl=quant.Salmon(index=testVars.testDir+"/salmonIndex/salIndex",transcriptome=testVars.cdna_big)
    assert sl.check_index()==True, "Failed salmon build_index"
    
    opts={"-o":testVars.testDir+"/salOut",
          "-l":"A",
          "-1":testVars.fq1,
          "-2":testVars.fq2}
    
    st=sl.run(subcommand='quant',target='tests/testout/salOut/quant.sf',**opts)
    assert st==True, "Failed to run salmon"
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    