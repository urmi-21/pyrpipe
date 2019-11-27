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

'''
#single end ERR2929684
#download sra->fq>qc
newSRA=sra.SRA('SRR5507343',testDir)
newSRA.downloadSRAFile()
newSRA.runFasterQDump(**{"-f":"","-t":testDir})

#run trimgalore
tg=qc.Trimgalore(**{"-j":"8","--length":"1"})  #specify to use 8 cores
bd=qc.BBmap()

#newSRA.performQC(bd)
newSRA.performQC(tg)

#run hisat and stie
hs=mapping.Hisat2(hisatInd)
hisatStatus=hs.runHisat2(newSRA,**{"-p":"10","--dta-cufflinks":""})
'''

'''
newSRA2=sra.SRA('ERR3527958',testDir)
newSRA2.downloadSRAFile()
newSRA2.runFasterQDump(**{"-f":""})
newSRA2.performQC(tg)
'''

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



#build pipeline with cleanup
#bbduk object
bd=qc.BBmap()
#trimgalore object to run trim_galore
tg=qc.Trimgalore(**{"-j":"8","--length":"1"})  #specify to use 8 cores
#Hisat2 object
hs=mapping.Hisat2(hisatInd)
#samtools object
samtOb=mapping.Samtools()
#Stringtie object
stieOb=assembly.Stringtie()

sraOb=sra.SRA('SRR5507495',testDir)

#download sra
sraOb.downloadSRAFile()
#run fastqdump;delete sra when done
sraOb.runFasterQDump(deleteSRA=True,**{"-f":"","-t":testDir})
#perform qc using trim_galore
sraOb.performQC(tg,deleteRawFastq=True)
#run hisat and store the status and return Sam file path
hisatStatus=hs.runHisat2(sraOb,**{"-p":"10","--dta-cufflinks":""})

#check if hisat is sucessful
if not hisatStatus[0]:
    raise Exception("ERROR: Hisat failed")
    
#remove qc corrected fastq
sraOb.deleteFastqFiles()
#run sam to sorted bam then run stringtie
stieOb.runStringtie(samtOb.samToSortedBam(hisatStatus[1],10,deleteSam=True,deleteOriginalBam=True),deleteInputBam=True,proc=10)






