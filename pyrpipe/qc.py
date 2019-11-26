#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 17:48:00 2019

@author: usingh
"""

from myutils import *

class Trimgalore:
    def __init__(self):
        self.programName="trim_galore"
        self.depList=[self.programName,'cutadapt']
        #check if hisat2 exists
        if not checkDep(self.depList):
            raise Exception("ERROR: "+ self.programName+" not found.")