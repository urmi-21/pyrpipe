#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 18:26:06 2020

@author: usingh
"""
from pyrpipe import quant
from pyrpipe import pyrpipe_utils as pu
from testingEnvironment import testSpecs
#import os

testVars=testSpecs()

def test_mapping():
    #create assemble object
    qob=quant.Quant()
    assert qob.category == "Quantification","Failed Quant init"
 
def test_kallisto():
    kl=quant.Kallisto(kallisto_index="")
    assert kl.check_index()==False, "Failed kallisto check_index"
    st=kl.build_index(index_path=testVars.testDir+"/kallistoIndex",index_name="kalIndex",fasta=testVars.cdna)
    assert st==True, "Failed to build kallisto index"
    opts={"-o":testVars.testDir+"/kalOut","--":(testVars.fq1,testVars.fq2)}
    st=kl.run_kallisto("quant",**opts)
    assert st==True, "Failed to run kallisto"
    
    
def test_salmon():
    sl=quant.Salmon(salmon_index="")
    assert sl.check_index()==False, "Failed salmon check_index"
    st=sl.build_index(index_path=testVars.testDir+"/salmonIndex",index_name="salIndex",fasta=testVars.cdna_big)
    assert st==True, "Failed to build salmon index"
    opts={"-o":testVars.testDir+"/salOut",
          "-l":"A",
          "-1":testVars.fq1,
          "-2":testVars.fq2}
    st=sl.run_salmon("quant",**opts)
    assert st==True, "Failed to run salmon"
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    