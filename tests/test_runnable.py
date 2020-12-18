#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 17:38:46 2020

@author: usingh
"""

from pyrpipe.runnable import Runnable
from testingEnvironment import testSpecs
import os

testVars=testSpecs()
cdna=testVars.cdna_big
workingDir=testVars.testDir
outdir=os.path.join(workingDir,'orfipy_out')
param={'--outdir':outdir,'--procs':'3','--dna':'orfs.fa','--bed':'orfs.bed'}


def test_basic():
    orfipy=Runnable(command='orfipy')
    target=os.path.join(outdir,'orfs.fa')
    orfipy.run(cdna,**param,target=target)
    
def test_targets():
    orfipy=Runnable(command='orfipy')
    targets=[os.path.join(outdir,'orfs.fa'),os.path.join(outdir,'orfs.bed')]
    orfipy.run(cdna,**param,target=targets)
    
def test_require():
    orfipy=Runnable(command='orfipy')
    req=cdna
    targets=[os.path.join(outdir,'orfs.fa'),os.path.join(outdir,'orfs.bed')]
    orfipy.run(cdna,**param,target=targets,requires=req)
    
def test_require2():
    orfipy=Runnable(command='orfipy')
    req=[cdna]
    targets=os.path.join(outdir,'orfs.bed')
    orfipy.run(cdna,**param,target=targets,requires=req)

