#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:21:01 2019

@author: usingh
"""

from myutils import *

class Stringtie:
    def __init__(self):
        self.programName="hisat2"
        #check if hisat2 exists
        if not checkDep([self.programName]):
            raise Exception("ERROR: "+ self.programName+" not found.")
            