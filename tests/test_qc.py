#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 28 01:18:25 2019

@author: usingh
"""

from pyrpipe import qc
from testingEnvironment import testSpecs
import os

def test_RNASeqQC():
    ob=qc.RNASeqQC()
    assert ob.category=="RNASeqQC", "RNASeqQC failed"