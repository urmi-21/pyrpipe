#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 19:06:34 2019

@author: usingh
"""

from pyrpipe import sra,mapping,assembly,qc,tools,pyrpipe_utils

workingDir="/home/usingh/work/urmi/hoap/test"

athalRunsSmol=['SRR1583780','SRR5507495','SRR5507442']
sraObjects=[]

for x in athalRunsSmol:
    thisSraOb=sra.SRA(x,workingDir)
    if thisSraOb.downloadSRAFile():
        sraObjects.append(thisSraOb)
    else:
        print("Download failed:"+x)

print("Following runs downloaded:")
for ob in sraObjects:
    print(ob.srrAccession)
    


for ob in sraObjects:
    ob.runFasterQDump(deleteSRA=True,**{"-e":"8","-f":"","-t":workingDir}) #use 8 threads

print("Fastq dump finished for:")
for ob in sraObjects:
    if ob.fastqFilesExistsLocally():
        print(ob.srrAccession)