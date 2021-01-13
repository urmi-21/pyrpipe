import sys
import pandas as pd
import os

input=sys.argv[1].split(',')
for f in input:
    print ('this f',f)

DIR=sys.argv[-1]

with open(input[0]) as f:
    thisdata=f.read().splitlines()
thisdata.pop(0)

names=[]
txids=[]
geneids=[]
for l in thisdata:
    thistx=l.split('\t')[0].split('|')[0]
    thisgene=l.split('\t')[0].split('|')[1]
    txids.append(thistx)
    geneids.append(thisgene)
df=pd.DataFrame({'TranscriptID':txids,'GeneID':geneids})

#read files in memory
for qf in input:
    name=qf.split('/')[-3]
    names.append(name)
    thisdata=pd.read_csv(qf,sep='\t',usecols=[3],skiprows=0)
    df[name]=thisdata['TPM']

#transcript counts
df_tx=df[['TranscriptID']+names].copy()
#write to file
df_tx.to_csv(DIR+'/results_TPM_tx.tsv',sep='\t',index=False)
#gene counts
df_gene=df[['GeneID']+names].copy()
df_gene = df_gene.groupby(['GeneID'],as_index = False).sum()
df_gene.to_csv(DIR+'/results_TPM_gene.tsv',sep='\t',index=False)
