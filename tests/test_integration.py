#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 10:27:44 2020

@author: usingh

Test various pyrpipe modules used with each other
"""

from pyrpipe import sra,qc,mapping,assembly,quant,tools
from pyrpipe import pyrpipe_utils as pu
from testingEnvironment import testSpecs
import os

testVars=testSpecs()
fq1=testVars.fq1
fq2=testVars.fq2
rRNAfasta=testVars.rRNAfa
srr='SRR978411'
workingDir=testVars.testDir
    
def test_pipeline1():
    sraOb=sra.SRA(srr,workingDir)
    st=sraOb.download_sra()
    assert st==True,"SRA download failed"
    
    st=sraOb.run_fasterqdump(delete_sra=False,**{"-e":"8","-f":"","-t":workingDir})
    assert st==True,"fqdump failed"
    
    bbdOpts={"ktrim":"r","k":"23","mink":"11","qtrim":"'rl'","trimq":"10","--":("-Xmx2g",),"ref":testVars.bbdukAdapters}
    bbdOb=qc.BBmap(**bbdOpts)
    st=sraOb.perform_qc(bbdOb)
    assert st==True,"bbduk failed"
    
    tgOpts={"--cores": "10", "-o":testVars.testDir, "--paired":"", "--": (fq1,fq2)}
    tg=qc.Trimgalore(**tgOpts)
    st=sraOb.perform_qc(tg)
    assert st==True,"tg failed"
    
    #runbowtie2
    bt=mapping.Bowtie2(bowtie2_index="")
    assert bt.check_index()==False, "Failed bowtie2 check_index"
    st=bt.build_index(testVars.testDir+"/btIndex","bowtieIndex",testVars.genome)
    assert st==True, "Failed to build bowtie2 index"
    st=bt.perform_alignment(sraOb)
    assert os.path.isfile(st)==True,"bowtie failed"
    
    hsOpts={"--dta-cufflinks":"","-p":"8"}
    hs=mapping.Hisat2(hisat2_index="",**hsOpts)
    st=hs.build_index(testVars.testDir,"hisatindex",testVars.genome)
    assert st==True, "Failed to build hisat2 index"
    #perform alignment with sraobject
    st=hs.perform_alignment(sraOb)
    assert os.path.isfile(st)==True,"hisat failed"
    
    hisatSam=st
    samOb=tools.Samtools(**{"-@":"8"})
    bam=samOb.sam_sorted_bam(hisatSam,delete_sam=False,delete_bam=False)
    assert os.path.isfile(bam)==True,"sam to bam failed"
    
    stie=assembly.Stringtie(reference_gtf=testVars.gtf)
    result=stie.perform_assembly(bam,out_dir=testVars.testDir)
    assert pu.check_files_exist(result)==True, "Failed stringtie"
    
    tr=assembly.Trinity()
    tr_out=tr.perform_assembly(sraOb,verbose=True)
    assert pu.check_files_exist(tr_out)==True, "Failed stringtie"

    kl=quant.Kallisto(kallisto_index="")
    assert kl.check_index()==False, "Failed kallisto check_index"
    st=kl.build_index(index_path=testVars.testDir+"/kallistoIndex",index_name="kalIndex",fasta=testVars.cdna)
    assert st==True, "Failed to build kallisto index"
    st=kl.perform_quant(sraOb)
    assert os.path.isdir(st)==True, "Failed to run kallisto"
    
    sl=quant.Salmon(salmon_index="")
    assert sl.check_index()==False, "Failed salmon check_index"
    st=sl.build_index(index_path=testVars.testDir+"/salmonIndex",index_name="salIndex",fasta=testVars.cdna)
    assert st==True, "Failed to build salmon index"
    
    st=sl.perform_quant(sraOb)
    assert os.path.isdir(st)==True, "Failed to run salmon"
    
    
    