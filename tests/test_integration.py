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

#srr='ERR3770564' #single end arabidopsis data
#srr='SRR978414' #small a thal paired end data
srr='SRR4113368'
workingDir=testVars.testDir
    
def test_pipeline1():
 
    sraOb=sra.SRA(srr,workingDir)
    st=sraOb.fastq_exists()
    assert st==True,"fasterq-dump failed"
    bbdOpts={"ktrim":"r","k":"23","mink":"11","qtrim":"'rl'","trimq":"10","ref":testVars.bbdukAdapters}
    bbdOb=qc.BBmap(None,**bbdOpts)
    st=sraOb.trim(bbdOb)
    assert st==True,"bbduk failed"
    
    
    tg=qc.Trimgalore()
    st=sraOb.qc(tg)
    assert st==True,"tg failed"
    
    #runbowtie2
    bt=mapping.Bowtie2(index=testVars.testDir+"/btIndex",genome=testVars.genome)
    assert bt.check_index()==True, "Failed bowtie2 check_index"
    st=bt.align(sraOb)
    assert os.path.isfile(st)==True,"bowtie failed"
    
    hsOpts={"--dta-cufflinks":"","-p":"8"}
    hs=mapping.Hisat2(index="")
    st=hs.build_index(testVars.testDir,"hisatindex",testVars.genome)
    assert st==True, "Failed to build hisat2 index"
    #perform alignment with sraobject
    st=hs.perform_alignment(sraOb,**hsOpts)
    assert os.path.isfile(st)==True,"hisat failed"
    
    hisatSam=st
    samOb=tools.Samtools(threads=2)
    bam=samOb.sam_sorted_bam(hisatSam,delete_sam=False,delete_bam=False)
    print(bam)
    assert os.path.isfile(bam)==True,"sam to bam failed"
    
    stie=assembly.Stringtie()
    result=stie.perform_assembly(bam,out_dir=testVars.testDir,reference_gtf=testVars.gtf)
    assert pu.check_files_exist(result)==True, "Failed stringtie"
    
    tr=assembly.Trinity()
    tr_out=tr.perform_assembly(sraOb,threads=2,verbose=True)
    assert pu.check_paths_exist(tr_out)==True, "Failed trinity"

    kl=quant.Kallisto(kallisto_index="")
    assert kl.check_index()==False, "Failed kallisto check_index"
    st=kl.build_index(index_path=testVars.testDir+"/kallistoIndex",index_name="kalIndex",fasta=testVars.cdna)
    assert st==True, "Failed to build kallisto index"
    st=kl.perform_quant(sraOb)
    print("asdsa"+st)
    assert os.path.isdir(st)==True, "Failed to run kallisto"
    
    sl=quant.Salmon(salmon_index="")
    assert sl.check_index()==False, "Failed salmon check_index"
    st=sl.build_index(index_path=testVars.testDir+"/salmonIndex",index_name="salIndex",fasta=testVars.cdna)
    assert st==True, "Failed to build salmon index"
    
    st=sl.perform_quant(sraOb,**{'--minAssignedFrags':'1'})
    assert os.path.isdir(st)==True, "Failed to run salmon"
    
    
    
