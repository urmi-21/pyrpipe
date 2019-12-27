#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 13:35:31 2019

@author: usingh
"""

class testSpecs:
    def __init__(self):
        self.srr='SRR1583780' #paired end Saccharomyces cerevisiae; RNA-Seq
        self.testDir="tests/testout"
        testfiles="tests/test_files/athaliana"
        self.fq1=testfiles+"/fastq/SRR978411_small_1.fastq"
        self.fq2=testfiles+"/fastq/SRR978411_small_1.fastq"
        self.genome=testfiles+"/genome/Arabidopsis_thaliana.TAIR10.dna.1and2.fa"
        self.gtf=testfiles+"/genome/Arabidopsis_thaliana.TAIR10.45_1and2.gtf"
        self.hisat2index=testfiles+"/hisat2index/athal"
        self.hisatSam=testfiles+"/mapping/hisat2.sam"
        self.hisatBam=testfiles+"/mapping/hisat2.bam"
        self.hisatSortedBam=testfiles+"/mapping/hisat2_sorted.bam"