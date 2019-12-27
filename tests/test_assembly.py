#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 15:56:20 2019

@author: usingh
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 16:24:44 2019

@author: usingh
"""
import pytest
from pyrpipe import assembly
from testingEnvironment import testSpecs
import os

testVars=testSpecs()

def test_assembly():
    srrID=testVars.srr
    testDir=testVars.testDir
    #create assemble object
    aob=assembly.Assembly()
    assert aob.category == "Assembler","Failed Assembly init"
    
def test_stringtie():
    input_bam=""