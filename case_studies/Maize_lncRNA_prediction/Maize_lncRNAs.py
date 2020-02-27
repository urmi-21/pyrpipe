#!/usr/bin/env python
# coding: utf-8

# In[3]:


from pyrpipe import sra,mapping,assembly,qc,tools
from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe
#First get the srr accessions of the runs. For this one can use the python package pysradb or R package sradb
#runs=['SRR3098746','SRR3098745','SRR3098744'] #from the study SRP068369
runs=['SRR765545'] #small test
#set up directories

workingDir="maize_out"
#create working directory
if not pu.check_paths_exist(workingDir):
    pu.mkdir(workingDir)


# ## Download Genome and GTF

# In[ ]:


GENOME=workingDir+"/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
GTF=workingDir+"/Zea_mays.B73_RefGen_v4.46.gtf"

if not pu.check_files_exist(GENOME):
    print("Downloading genome fasta file")
    wget="wget ftp://ftp.ensemblgenomes.org/pub/release-46/plants/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz -q -O "+GENOME+".gz"
    pe.execute_command(wget.split(),verbose=True,logs=False)
    pe.execute_command(['gunzip',GENOME+".gz"],verbose=True,logs=False)
    
if not pu.check_files_exist(GTF):
    print("Downloading GTF file")
    wget="wget ftp://ftp.ensemblgenomes.org/pub/release-46/plants/gtf/zea_mays/Zea_mays.B73_RefGen_v4.46.gtf.gz -q -O "+GTF+".gz"
    pe.execute_command(wget.split(),verbose=True,logs=False)
    pe.execute_command(['gunzip',GTF+".gz"],verbose=True,logs=False)


# ## Download data, pre-process

# In[5]:


sraObjects=[]
for x in runs:
    thisSraOb=sra.SRA(x,workingDir)
    if thisSraOb.download_sra():
        sraObjects.append(thisSraOb)
    else:
        print("Download failed:"+x)
        
#perform fastq dump and qc

#create a Trimgalore object
tgOpts={"--cores": "10"}
tg=qc.Trimgalore(**tgOpts)
#NOTE: To download fastq directly, instaead of .sra, one can use the download_fastq() method
for x in sraObjects:
    #to fastq
    x.run_fasterqdump(delete_sra=True,**{"-e":"20","-f":"","-t":workingDir}) #use 20 threads
    #perform qc using trim galore
    x.perform_qc(tg)
    

        


# ## Map using STAR

# In[ ]:


starParams={"--outFilterType":"BySJout",
            "--runThreadN":"8",
            "--outSAMtype": "BAM SortedByCoordinate"
            }

star=mapping.Star(star_index="",**starParams) #provided index is invalid

#create star index
indexOut=workingDir+"/starindex"
inFasta=GENOME
star.build_index(indexOut,inFasta)


# ## Transcript assembly using StringTie

# In[ ]:


#Create object for stringtie. This will be used for all the bam files.
st=assembly.Stringtie(reference_gtf=GTF)
gtfList=[]
for x in sraObjects:
    star_out_dir=star.perform_alignment(x,objectid=x.srr_accession)
    bam=star_out_dir+"/Aligned.sortedByCoord.out.bam"
    gtfList.append(st.perform_assembly(bam,objectid=x.srr_accession,**{"-p":"25"}))   

print(gtfList)


# ## lncRNA prediction using PLncPRO
# We will use [PLncPRO](https://github.com/urmi-21/PLncPRO) for prediction of lncRNAs. Currently, PLncPRO is not integrated into `pyrpipe` so we will use the `pyrpipe_engine` module directly to execute.

# In[ ]:


#import pyrpipe modules
from pyrpipe import pyrpipe_engine as pe
#install plncpro
pe.execute_command("pip install plncpro".split(),verbose=True,quiet=False,logs=False)
#OR
#!pip install plncpro


genome="maize_data/Zea_mays.B73_RefGen_v4.dna.toplevel.1_10.fa"
model="monocot_model/monocot.model"
blastdb="uniprot/uniprotdb"
for i in range(len(gtfList)):
    thisOb=sraObjects[i]
    #first extract transcripts using gffread
    tx_file=thisOb.location+"/transcripts.fa"
    cmd="gffread -w "+tx_file+" -g maize_data/Zea_mays.B73_RefGen_v4.dna.toplevel.1_10.fa "+gtfList[i]
    pe.execute_command(cmd.split(" "),verbose=False,quiet=False,logs=True,objectid=thisOb.srr_accession,command_name="gffread")
    
    #Optional step use biopython to filter transcripts by len
    #out_file=thisOb.location+"/transcripts_filter.fa"
    #output_handle = open(out_file, "w")
    #for record in SeqIO.parse(tx_file, "fasta"):
        # keep tx between 200 and 1000
    #    if len(record)>=500 and len(record)<=1000:
    #        #write to temp file
    #        SeqIO.write(record, output_handle, "fasta")

    
    #run plncpro
    outdir=thisOb.location+"/plncpro_out"
    outfile="plncpro_predictions"
    cmd="plncpro predict -i "+tx_file+" -o "+outdir+" -p "+outfile+" -t 25 --min_len 200 -d "+blastdb+" -m "+model+" -v -r"
    pe.execute_command(cmd.split(),verbose=False,quiet=False,logs=True,objectid=thisOb.srr_accession,command_name="plncpro predict")
        


# ## Generate reports

# In[ ]:


#NOTE: Following commands are executed in shell, hence the ! before each command
#get_ipython().system('pyrpipe_diagnostic.py report pyrpipe_logs/2020-01-22-18_14_47_pyrpipe.log')
#get_ipython().system('pyrpipe_diagnostic.py benchmark pyrpipe_logs/2020-01-22-18_14_47_pyrpipe.log')
#get_ipython().system('pyrpipe_diagnostic.py shell pyrpipe_logs/2020-01-22-18_14_47_pyrpipe.log')
#get_ipython().system('pyrpipe_diagnostic.py multiqc -o ./multiqc_report pyrpipe_logs/2020-01-22-18_14_47_pyrpipe.log')


# In[ ]:




