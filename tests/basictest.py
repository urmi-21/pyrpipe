#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:53:26 2019

@author: usingh
"""

from pyrpipe import sra,mapping,assembly,qc,tools,pyrpipe_session


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

sraOb=sra.SRA('SRR1583780',testDir)

#download sra
sraOb.downloadSRAFile()
#run fastqdump;delete sra when done
sraOb.runFasterQDump(deleteSRA=True,**{"-f":"","-t":testDir})
#perform qc using trim_galore
sraOb.performQC(tg,deleteRawFastq=True)
#run hisat and store the status and return Sam file path
hisatSam=hs.runHisat2(sraOb,**{"-p":"10","--dta-cufflinks":""})

#check if hisat is sucessful
if not hisatSam:
    raise Exception("ERROR: Hisat failed")
    
#remove qc corrected fastq
sraOb.deleteFastqFiles()
#run sam to sorted bam then run stringtie
gtfS=stieOb.runStringtie(samtOb.samToSortedBam(hisatSam,10,deleteSam=True,deleteOriginalBam=True),deleteInputBam=True,proc=10)
"""



btIndex="/home/usingh/work/urmi/hoap/test/bowtieIndex/rRNAindex"
#riboseq SRR3590744
sraOb=sra.SRA('SRR5507495',testDir)
#download sra
sraOb.downloadSRAFile()
#run fastqdump;delete sra when done
sraOb.runFasterQDump(deleteSRA=False,**{"-f":"","-t":testDir})

tgOb=qc.Trimgalore()

#sraOb.performFastqQC(tgOb)
pathToAdapters="/home/usingh/lib_urmi/softwares/bbmap/resources/adapters2.fa"
bbdOpts={"ktrim":"r","k":"23","mink":"11","qtrim":"'rl'","trimq":"10","--":("-Xmx2g",),"ref":pathToAdapters}
bbdOb=qc.BBmap(**bbdOpts)

#sraOb.performFastqQC(bbdOb)
#status=bbdOb.performCleaning(sraOb,"/home/usingh/work/urmi/hoap/test/bowtieIndex/euk_combined_rRNA.fa")

#print(status)


#run bbmap
#bd=qc.BBmap()
#sraOb.performQC(bd)

#run bowtie
#bob=mapping.Bowtie2(btIndex)
#unMappedReads=bob.runBowTie2(sraOb)

#update fastq as
#sraOb.localfastqPath=unMappedReads


#build hisat index

hsOpts={"--dta-cufflinks":"","-p":"12","--mp": "1,1", "--no-spliced-alignment":"", "--rdg": "10000,10000", "--rfg": "10000,10000"}
hs=mapping.Hisat2(hisat2Index="/home/usingh/work/urmi/hoap/test/yeastInd2/index22",**hsOpts)
#hsbArgs={"-p":"8","-a":"","-q":""}
#if hs.buildHisat2Index("/home/usingh/work/urmi/hoap/test/yeastInd2","index22","/home/usingh/work/urmi/hoap/test/hisatYeast/S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113.fsa",**hsbArgs):
#    print("Success")
    
#run hisat
sam=hs.performAlignment(sraOb,**{"--dta-cufflinks":"","-p":"8"})

#get sorted bam
samOb=tools.Samtools(**{"-@":"8"})
bam=samOb.samToSortedBam(sam,deleteSam=True,deleteOriginalBam=True)


#bt2=mapping.Bowtie2("/home/usingh/work/urmi/hoap/test/bowtieIndex/rRNAindex")
#bt2.performAlignment(sraOb)

#run stringtie
#st=assembly.Stringtie()
#g1=st.performAssembly(bam)


#gtfs=(g1,)
#test stmerge
#merged=st.performStringtieMerge(g1,g1,outFileSuffix="_stOUT",overwrite=True)
#if not merged:
#    print("Fail")
"""
bamList=[]
sraObList=[]
for s in ['SRR1583780','SRR5507495','SRR5507442','SRR5507362']:
    sraOb=sra.SRA(s,testDir)
    #download sra
    sraOb.downloadSRAFile()
    #run fastqdump;delete sra when done
    sraOb.runFasterQDump(deleteSRA=True,**{"-f":"","-t":testDir})
    sraObList.append(sraOb)
    sam=hs.performAlignment(sraOb,**{"--dta-cufflinks":"","-p":"8"})
    #get sorted bam
    bam=samOb.samToSortedBam(sam,deleteSam=True,deleteOriginalBam=True)
    bamList.append(bam)
""" 

#bam merge
#mergedBam=samOb.mergeBamFiles(*bamList,outPath=testDir,outFileName="myMergedXXDD",*{"-f":""})

#using pysam
#import pysam
#pysam.merge("-@","8","myMerge",*bamList,"-f")


#portc
#pob=tools.Portcullis()
#refGenome="/home/usingh/work/urmi/hoap/test/hisatYeast/S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113.fsa"
#portDir=pob.runPortcullisFull(refGenome,mergedBam,outDir=testDir+"/portOut",deleteOriginalBamFile=True)


#save work space
#pyrpipe_session.savePyrpipeWorkspace(filename="sess",outDir=testDir)

#test star
"""
starParams={"--outFilterType":"BySJout",
            "--runThreadN":"8",
            "--outSAMtype": "BAM SortedByCoordinate"
            }
star=mapping.Star(starIndex="/home/usingh/work/urmi/hoap/test/yeaststarIndex",**starParams)

starOut=testDir+"/starout"

star.performAlignment(sraOb)

star2=mapping.Star(**starParams)
starIndOut="/home/usingh/work/urmi/hoap/test/si2"
inFasta="/home/usingh/work/urmi/hoap/test/hisatYeast/S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113.fsa"
sind=star2.buildStarIndex(starIndOut,inFasta,**{"--genomeSAindexNbases":"8","--outFileNamePrefix":testDir})

print (sind)

#run cufflinks
refGTF="/home/usingh/work/urmi/hoap/test/hisatYeast/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.gff"
cl=assembly.Cufflinks()

clout=cl.performAssembly(bam,**{"--num-threads":"28","--no-update-check":""})
print(clout)
"""
#runribocode
#rbc=tools.Ribocode()
#rbc.runRibocode()




#test mikado
atRef="/home/usingh/work/urmi/hoap/test/athalData/ref/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
mk=tools.Mikado()
#find and save gtf list
#mklist=mk.searchGTFtolist("athalGtfList",searchPath="/home/usingh/work/urmi/hoap/test/athalData/sraData",outDir="/home/usingh/work/urmi/hoap/test/athalData/sraData")

mikadoDir="/home/usingh/work/urmi/hoap/test/mikadoTutorial"
junc=mikadoDir+"/junctions.bed"
ref=mikadoDir+"/chr5.fas"
#Make sure the paths in list file are global.
mklist=mikadoDir+"/list.txt"
scoring=mikadoDir+"/plants.yaml"
config=mk.runMikadoConfigure(mklist,ref,"permissive",scoring,junc,"mkConfig")
print(config)

#run mikado
mk.runMikadoPrepare(config,outDir="/home/usingh/work/urmi/hoap/test/sample_data/mikadoPrepOut")

mk.runMikadoSerialise(config,mikadoDir+"/mikadoPrepOut/mikado_prepared.fasta"+"",junc,mikadoDir+"/uniprot_sprot_plants.fasta",mikadoDir+"/mikado.bed",mikadoDir+"/mikado.blast.xml.gz",outDir=mikadoDir+"/serOut")


#bamList=searchFilesLocally("/home/usingh/work/urmi/hoap/test/athalData/sraData","*.bam")
#bamList=["/home/usingh/work/urmi/hoap/test/athalData/sraData/SRR971778/SRR971778_hisat2_sorted.bam",
         #"/home/usingh/work/urmi/hoap/test/athalData/sraData/SRR978411/SRR978411_hisat2_sorted.bam",
         #"/home/usingh/work/urmi/hoap/test/athalData/sraData/SRR976159/SRR976159_hisat2_sorted.bam"
         #]

#mergedBam=samOb.mergeBamFiles(*bamList,outPath="/home/usingh/work/urmi/hoap/test/athalData/sraData/",outFileName="atMerged",*{"-f":""})
#run portculis
#pob=tools.Portcullis()
#portOut=pob.runPortcullisFull(atRef,mergedBam,outDir=testDir+"/home/usingh/work/urmi/hoap/test/athalData/sraData/portOut")

#run configure

#atJunc=""
#

