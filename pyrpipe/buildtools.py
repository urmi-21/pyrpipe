#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 16:12:47 2021

@author: usingh
"""

import os


def build_tools():
    cmd_ins_mamba="conda install -c conda-forge mamba"
    #cmd_update_py='mamba install python=3.8'
    cmd_bioconda='mamba install -c bioconda sra-tools=2.10.9 trim-galore=0.6.6  bbmap=38.87 hisat2=2.2.1  star=2.7.7a bowtie2=2.4.2 kallisto=0.46.2 salmon stringtie=2.1.4 orfipy=0.0.3 samtools=1.11'
    
    print(cmd_ins_mamba)
    os.system(cmd_ins_mamba)
    print(cmd_bioconda)
    os.system(cmd_bioconda)
    
    
    
    
    