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
        self.fq2=testfiles+"/fastq/SRR978411_small_2.fastq"
        self.genome=testfiles+"/genome/Arabidopsis_thaliana.TAIR10.dna.1and2.fa"
        self.gtf=testfiles+"/genome/Arabidopsis_thaliana.TAIR10.45_1and2.gtf"
        self.cdna=testfiles+"/genome/Arabidopsis_thaliana.TAIR10.cdna.5000.fa"
        self.hisat2index=testfiles+"/hisat2index/athal"
        self.hisatSam=testfiles+"/mapping/hisat2.sam"
        self.hisatBam=testfiles+"/mapping/hisat2.bam"
        self.hisatSortedBam=testfiles+"/mapping/hisat2_sorted.bam"
        self.rRNAfa="tests/test_files/euk_1000_rRNA.fa"
        self.starSortedBam=testfiles+"/mapping/Aligned.sortedByCoord.out.bam"
        
        self.mikadofiles="tests/test_files/athaliana/mikado"
        
        self.bbdukAdapters="tests/test_files/adapters2.fa"
        self.uniprot="tests/test_files/uniprot_sprot_plants_1000.fasta"
        
        