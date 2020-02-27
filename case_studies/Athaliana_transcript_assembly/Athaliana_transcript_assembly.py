#!/usr/bin/env python
# coding: utf-8

# # Analysis of *A. thaliana* RNA-Seq data with pyrpipe
# Use A thaliana public RNA-Seq data to assemble transcripts.

# In[1]:


from pyrpipe import sra,mapping,assembly,qc,tools
from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe
#First get the srr accessions of the runs. For this one can use the python package pysradb or R package sradb
#i will consider following randomly selected accessions
#athalRuns=['SRR976159','SRR978411','SRR978410','SRR971778','SRR1058116','SRR1058118','SRR1058121','SRR1058110','SRR1058120','SRR1058117','SRR1104134','SRR1104133','SRR1104135','SRR1104136','SRR1105825']
athalRunsSmol=['SRR976159','SRR978411','SRR971778']
#set your working directory if you don't want to use the current working directory
workingDir="athal_out"
#create working directory
if not pu.check_paths_exist(workingDir):
    pu.mkdir(workingDir)


# ## Download genome and gtf

# In[ ]:


GENOME=workingDir+"/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
GTF=workingDir+"/Arabidopsis_thaliana.TAIR10.45.gtf"

if not pu.check_files_exist(GENOME):
    print("Downloading genome fasta file")
    wget="wget ftp://ftp.ensemblgenomes.org/pub/release-46/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz -q -O "+GENOME+".gz"      
    pe.execute_command(wget.split(),verbose=True,logs=False)
    pe.execute_command(['gunzip',GENOME+".gz"],verbose=True,logs=False)

if not pu.check_files_exist(GTF):
    print("Downloading GTF file")
    wget="wget ftp://ftp.ensemblgenomes.org/pub/release-46/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.46.gtf.gz -O "+GTF+".gz"
    pe.execute_command(wget.split(),verbose=True,logs=False)
    pe.execute_command(['gunzip',GTF+".gz"],verbose=True,logs=False)




# ## Download data and create SRA objects
# First we can download all data to disk and create pyrpipe.SRA objects

# In[3]:



##download all data in athalRuns
sraObjects=[]

for x in athalRunsSmol:
    thisSraOb=sra.SRA(x,workingDir)
    if thisSraOb.download_sra():
        sraObjects.append(thisSraOb)
    else:
        print("Download failed:"+x)

#NOTE: To download fastq directly, instaead of .sra, one can use the download_fastq() method
print("Following runs downloaded:")
for ob in sraObjects:
    print(ob.srr_accession)


# ## Saving current session
# A reason why I have first downloaded the SRA files is that **in a typical HPC setting, one might have access to special data-transfer nodes**. These nodes could be used for downloading data efficiently but does not allow expensive computations. On the other hand data could also be downloaded from compute nodes **but you will burn most of your computing time/allocations for only downloading the data**. Thus it might be a good idea to download data separately and then start the processing.    #
# We can save the objects created with pyrpipe and restore our session later on a compute node.

# In[ ]:


# save current session
#from pyrpipe import pyrpipe_session
#pyrpipe_session.save_session(filename="mySession",add_timestamp=True,out_dir=workingDir)


# ## Restoring saved session
# We can restore the pyrpipe session using the saved session file (saved with .pyrpipe extension).
#
# **Note** After restoring session a new log file will generated to store the logs.

# In[ ]:


#first clear current session used by notebook
#get_ipython().run_line_magic('reset', '')
#print(sraObjects)


# In[ ]:


#restore session
#from pyrpipe import pyrpipe_session
#update the pyrpipe session file below
#pyrpipe_session.restore_session("athal_out/mySession_20200129112604.pyrpipe")
#print(sraObjects)


# ## Processing sra files
#
#  After restoring session we can proceed. So far we have downloaded data the sra files and sraObjects contsins all the SRA objects coressponding to each SRR Accession.
#

# ## Convert sra to fastq file
# We can convert .sra files to .fastq. Using ```delete_sra=True``` will delete the downloaded sra file from disk. Note the second argument ```**{"-e":"8","-f":"","-t":workingDir}``` is basically a dict containing additional fasterq-dump parameters.

# In[4]:


for ob in sraObjects:
    ob.run_fasterqdump(delete_sra=True,**{"-e":"8","-f":"","-t":workingDir}) #use 8 threads

print("Fastq dump finished for:")
for ob in sraObjects:
    if ob.fastqFilesExistsLocally():
        print(ob.srr_accession)


# ## Performing fastq quality control
# After running fasterq-dump, the fastq files will be updated in each SRA object. To perform fastq quality control we can use ```trimgalore``` or ```bbduk.sh```.     
# In[5]:


#using bbduk
pathToAdapters="adapters2.fa"
#arguments to pass to bbduk
bbdOpts={"ktrim":"r","k":"23","mink":"11","qtrim":"'rl'","trimq":"10","--":("-Xmx2g",),"ref":pathToAdapters}
#an object for running bbduk.sh with specified parameters
bbdOb=qc.BBmap(**bbdOpts)
#start QC
for ob in sraObjects:
    ob.perform_qc(bbdOb)

#after finishing view the current fastq files in the sra objects

for ob in sraObjects:
    print("SRR Accession: {}, fastq files: {}. {}".format(ob.srr_accession,ob.localfastq1Path,ob.localfastq2Path))

    if ob.fastqFilesExistsLocally():
          print("Both files exist!!")
    else:
          print("Error")
          raise Exception("Fastq files not found")


# ## Aligning clean reads to the reference genome
# After finishing fastq quality control we will map reads to the reference genome.

# In[6]:


#using hisat2
hsOpts={"--dta-cufflinks":"","-p":"8"}
hs=mapping.Hisat2(hisat2_index="",**hsOpts)


# In[7]:


#We can build hisat2 index if one doesnt already exist. This index will be bound to the Hisat2 object, hs.
hisat2_buildArgs={"-p":"8","-a":"","-q":""}
#start building
#parameters are out directory, index name, reference genome
if hs.build_index(workingDir+"/athalIndex","athalInd",GENOME,**hisat2_buildArgs) :
    print("Indexing done.")

#check the index present in hisat2 object
if hs.check_index():
    print("Index {} exists".format(hs.hisat2_index))



# In[8]:


#start alignment
samList=[]
for ob in sraObjects:
    print("Processing {}...".format(ob.srr_accession))
    thisSam=hs.perform_alignment(ob,**{"-p":"10"}) #note parametrs supplied here will replace existing parameters passed during object construction
    if thisSam:
        samList.append(thisSam)
print("Alignment done!! Sam files:"+ ",".join(samList))


# ## Using samtools
# ```pyrpipe``` implemnts a basic high-level samtools API through which samtools functionality could be accessed. Note that users can also use the library ```pysam``` to get advance SAM/BAM/VCF/BCF functionality.

# In[9]:


samOb=tools.Samtools(**{"-@":"8"})
#sam to sorted bam
bamList=[]
i=0
for sam in samList:
    print("Processing:"+sam)
    thisBam=samOb.sam_sorted_bam(sam,delete_sam=True,delete_bam=True,objectid=sraObjects[i].srr_accession) #add the object id to keep track of process and object. helpful in debugging and reports later
    i+=1
    if thisBam:
        bamList.append(thisBam)
print("Sorted bam files:"+",".join(bamList))

###Some Examples using pysam###
#for details see: https://pysam.readthedocs.io/en/latest/
#import pysam
#pysam.sort("-@","8","-o","sortedBam.bam","in.bam)
#pysam.merge("-@","8","myMerge",*bamList,"-f")


# ## Transcript assembly
# We can use stringtie to perform transcript assembly.

# In[10]:


st=assembly.Stringtie(reference_gtf=GTF)
gtfList=[]
i=0
for bam in bamList:
    print("Processing:"+bam)
    gtfList.append(st.perform_assembly(bam,objectid=sraObjects[i].srr_accession))
    i+=1

print("Final GTFs:"+",".join(gtfList))


# ## Generating analysis reports
# pyrpipe_diagnostic.py lets user generate different types of reports and summaries. Following commands can be run from shell.
#
#
# **Generate a pdf report**
# ```$ pyrpipe_diagnostic.py report pyrpipe_logs/2019-12-24-16_14_55_pyrpipe.log```
#
# [Output](https://github.com/urmi-21/pyrpipe/blob/master/examples%28case-studies%29/Athaliana_transcript_assembly/2020-01-29-12_00_45_pyrpipe.pdf)
#
# ***Dump all commands to a shell file***
# ```pyrpipe_diagnostic.py shell pyrpipe_logs/2019-12-24-16_14_55_pyrpipe.log```
#
# [Output](https://github.com/urmi-21/pyrpipe/blob/master/examples%28case-studies%29/Athaliana_transcript_assembly/2020-01-29-12_00_45_pyrpipe.sh)
#
#
# **Generate multiqc report**
# ```$ pyrpipe_diagnostic.py multiqc -r pyrpipe_logs/2019-12-24-16_14_55_pyrpipe.log```
#
# [Output](https://github.com/urmi-21/pyrpipe/blob/master/examples%28case-studies%29/Athaliana_transcript_assembly/multiqc_report.html)
#
#
# **Generate runtime benchmarks**
# ```pyrpipe_diagnostic.py benchmark pyrpipe_logs/2019-12-24-16_14_55_pyrpipe.log```
#
# [Output](https://github.com/urmi-21/pyrpipe/tree/master/examples%28case-studies%29/Athaliana_transcript_assembly/benchmark_reports)
#

# In[ ]:

