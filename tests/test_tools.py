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
junctions=""


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


#TODO check portcullis installation in travis    
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
 

def test_mikado_full():
    mk=tools.Mikado()
    gtfdir=testVars.mikadofiles
    out_dir=testVars.testDir+"/mikadoout"
    junc=testVars.testDir+"/3-filt/portcullis_filtered.pass.junctions.bed"
    list_file=mk.createMikadoGTFlist("mikadolist",out_dir,gtfdir)
    #diamond obj to use
    dm=tools.Diamond(index="",mode='sensitive')
    dm.build_index(testVars.uniprot,"diamondDB",out_dir=testVars.testDir+"/dout")
    mode="permissive"
    scoring="plants.yaml"
    mkout=mk.runMikadoFull(list_file,testVars.genome,mode,scoring,junc,"mkconf",testVars.uniprot,dm,out_dir=testVars.testDir+"/mikadoout",verbose=False)
    st=pu.check_paths_exist(mkout)
    assert st==True, "Mikado Full failed" 

 
def test_diamond():
    dm=tools.Diamond(index="",mode='sensitive')
    dm.build_index(testVars.uniprot,"diamondDB",out_dir=testVars.testDir+"/dout")
    dm.run_align(testVars.cdna, "diamond_out", command="blastx", out_fmt=6, fmt_string="qseqid sseqid evalue pident", out_dir=testVars.testDir+"/dout")
    st=pu.check_files_exist(testVars.testDir+"/dout/diamond_out")
    assert st==True, "Diamond failed" 
    

def test_transdecoder():
    td=tools.Transdecoder()
    longOrfOut=td.run_transdecoder_longorfs(testVars.bbdukAdapters,out_dir=testVars.testDir+"/longorfsout")
    preddir=testVars.testDir+"/predout"
    predout=td.run_transdecoder_predict(testVars.bbdukAdapters,longOrfOut,out_dir=preddir)
    st=pu.check_paths_exist(predout)
    assert st==True, "TransDecoder failed" 
    

    
    
    
    
    