#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:53:26 2019

@author: usingh
"""

import sra,mapping,assembly,qc


#########define directories indices, reference gtf etc####
testDir="/home/usingh/work/urmi/hoap/test"
hisatInd="/home/usingh/work/urmi/hoap/test/hisatYeast/S288C_reference_genome_R64-2-1_20150113/yeastIndex"


#download sra->fq>qc
newSRA=sra.SRA('SRR1583780',testDir)
newSRA.downloadSRAFile()
newSRA.runFasterQDump()

#run trimgalore
tg=qc.Trimgalore()
tg.runTrimGalorePaired()


"""
yeastList=['SRR1583780','SRR5507495','SRR5507442','SRR5507362','SRR5507343','SRR5507356','SRR5507413','SRR5507339','SRR5507399','SRR5507353','SRR5507415','SRR5507444','SRR5507419','SRR5507379','SRR5507434']
    

#hisat object
hs=mapping.Hisat2(hisatInd)
#samtools object
samtOb=mapping.Samtools()
#string tie object
stieOb=assembly.Stringtie()

#trim_galore object
tg=qc.Trimgalore()

#bbmap object
bm=qc.BBmap()

    
    
shortList=yeastList[0:10]
sraOBs=[]
#download all sra
for r in shortList:
    thisOb=sra.SRA(r,testDir)
    thisOb.downloadSRAFile()
    thisOb.runFasterQDump(**{"-e":"8","-S":"","--skip-technical":"","-t":testDir,"-f":""})
    #run hisat
    hisatStatus=hs.runHisat2(thisOb,**{"-p":"10","--dta-cufflinks":""})
    if hisatStatus[0]:
        temp=samtOb.samToSortedBam(hisatStatus[1],10)
        stieOb.runStringtie(temp,8)
        
    else:
        print("Failed Hisat for "+r)
""" 