#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 16:24:44 2019

@author: usingh
"""

from pyrpipe import sra


print ("testing sra...")
newOb=SRA('SRR10408795')
print(newOb.getSrrAccession())
print(newOb.sraFileExistsLocally())
newOb.downloadSRAFile(**{"-O": "/home/usingh/work/urmi/hoap/test", "Attr2": "Val2","-q":""})
print(newOb.sraFileExistsLocally())

print("done")