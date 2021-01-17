#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 15:24:42 2021

@author: usingh
"""

from pyrpipe import sra
from pyrpipe import pyrpipe_utils as pu
import os


def runtest():
    failed=False
    sraob=sra.SRA('ERR726985',directory='./pyrpipe_sratest')
    if not sraob.fastq_exists(): pu.print_boldred('Test failed'); failed=True
    pu.print_notification('Cleaning up...')
    sraob.delete_fastq()
    os.rmdir(sraob.directory)
    
    if failed:
        pu.print_boldred('Paired end test failed'); 
        failed=False
    
    sraob=sra.SRA('SRR2134545',directory='./pyrpipe_sratest')
    if not sraob.fastq_exists(): pu.print_boldred('Test failed'); failed=True
    pu.print_notification('Cleaning up...')
    sraob.delete_fastq()
    os.rmdir(sraob.directory)
    
    if failed:
        pu.print_boldred('Single end test failed'); 
        failed=False
        
    if not failed:
        pu.print_green('\n#####################All Tests Passed#####################\n')
        os.rmdir('./pyrpipe_sratest')