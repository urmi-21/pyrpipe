#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GTEx processing pipeline
@author: usingh

python processGTEX.py <infile> <analysis_type>
in file contains gtex sample ids e.g.
GTEX-1117F-0226-SM-5GZZ7
GTEX-1117F-0426-SM-5EGHI
GTEX-1117F-0526-SM-5EGHJ

Example:
#to align and transcript assembly
python processGTEX.py gtex_ids align

#to run quant
python processGTEX.py gtex_ids quant

"""

import sys,os,yaml
from pyrpipe import pyrpipe_engine as pe
from pyrpipe import pyrpipe_utils as pu
from pyrpipe import sra,mapping,assembly,tools,quant
from pyrpipe.runnable import Runnable
from pyrpipe import _threads,_force,_dryrun



class Bamtofastq(Runnable):
    def __init__(self,*args,threads=None,**kwargs):
        super().__init__(*args,command='bamtofastq',**kwargs)
        self._deps=[self._command]
        self._param_yaml='bamtofastq.yaml'
        self._args_style='JAVA'
        self.resolve_parameter("--procs",threads,_threads,'_threads')

    def bamtofq(self,bam,oid,rm_bam=True):
        out_dir=pu.get_file_directory(bam)
        fastq_name=os.path.join(out_dir,oid)
        sname=os.path.join(out_dir,'s.fq')
        oname=os.path.join(out_dir,'o.fq')
        o2name=os.path.join(out_dir,'o2.fq')
        #tempfilename
        tmpdir=os.environ.get('LOCAL')
        if not tmpdir: tmpdir='./'
        tmpfile=os.path.join(tmpdir,pu.get_file_basename(bam)+pu.get_timestamp())

        internal_kwargs={'F':fastq_name+'_1.fastq',
                         'F2':fastq_name+'_2.fastq',
                         'S':sname,
                         'O':oname,
                         'O2':o2name,
                         'T':tmpfile,
                         'filename':bam}

        #call run
        status=self.run(None,objectid=oid,**internal_kwargs)
        if status and rm_bam: pe.delete_file(bam)
        return status

#####################Functions#############################################
@pe.dryable
def check_paired(bam):
    header=pe.get_shell_output(['samtools', 'view',bam,'-H'],verbose=False)
    if not header[0] == 0:
        print("Invalid BAM file")
        return False
    header=header[1]
    for l in header.split('\n'):
        if "@PG" in l:
            temp=l.split()
            index=temp.index('--readFilesIn')
            if (not temp[index+1].startswith('--')) and (not temp[index+2].startswith('--')):
                return True
            return False
    return False

def sortbam(bam,oid):
        outfile=pu.get_file_basename(bam)+"_sorted.bam"
        outdir=pu.get_file_directory(bam)
        outpath=os.path.join(outdir,outfile)
        cmd='sambamba sort -t 25 -m 100G -o '+outpath+' '+ bam 
        st=pe.execute_command(cmd.split(),logs=True,objectid=oid)
        if not st: 
                return ""
        return outpath


#####################################################################################

#argv[1]: file containg gtex/tcga ids
idsfile=sys.argv[1]
analysis=sys.argv[2]
runquant=False
runalign=False
if analysis=='quant': runquant=True
if analysis=='align': runalign=True

with open(idsfile) as f:
    data=f.read().splitlines()
#set infile dir as workdir
basedir=pu.get_file_directory(idsfile)


#pyrpipe objects
star=mapping.Star()
#Create stringtie object
stieobj=assembly.Stringtie()
#biobambam
biobb=Bamtofastq()
#salmon for quant
salmon=quant.Salmon()

#delete final sorted bam 
delete_bam=True

#out_dir is same name as input file
out_dir=basedir
progresslog=open(idsfile+'_progress.log','w')
bam_suffix='.Aligned.sortedByCoord.out.patched.md.bam'
#process each bam
for line in data:
    temp=line.split('\t')
    sampleid=temp[0]
    #add fid to out dir
    this_outdir=os.path.join(basedir,sampleid)
    pu.mkdir(this_outdir)
    bam=os.path.join(this_outdir,sampleid+bam_suffix)
    bam_dir=pu.get_file_directory(bam)

    #print('bam path:'+bam)
    if not pu.check_files_exist(bam):
        print(bam+' not found.')
        error=sampleid+'\tbam file not found'
        progresslog.write(error+'\n')
        continue
    #print('check paired') 
    if not check_paired(bam):
        print(bam+' is not paired. Deleting...')
        pe.delete_file(bam)
        error=sampleid+'\tbam file not paired'
        progresslog.write(error+'\n')
        continue
    
    #convert bam to fastq
    st=biobb.bamtofq(bam,sampleid)
    if not st:
        print('bam to fastq failed.')
        #log this message
        error=sampleid+'\tbam to fastq failed'
        progresslog.write(error+'\n')
        continue

    #create SRA object
    fq1=os.path.join(bam_dir,sampleid+'_1.fastq')
    fq2=os.path.join(bam_dir,sampleid+'_2.fastq')
    sobj=sra.SRA(fastq=fq1,fastq2=fq2)
    #check files exist
    if not _dryrun and not sobj.fastq_exists():
        print('Error creating SRA obj')
        error=sampleid+'\tFastq not found'
        progresslog.write(error+'\n')
        continue
    
    if runalign:
        #perform alignment
        #start star alignment
        star_bam=starobj.perform_alignment(sobj) #returns bam path
        #delete fastq
        sobj.delete_fastq()
        if not star_bam:
            print('Star failed for '+bam)
            #log
            error=sampleid+'\tSTAR alignment failed'
            progresslog.write(error+'\n')
            continue

        sorted_bam=sortbam(star_bam,sampleid)
        #delete original bam
        pe.delete_file(star_bam)
        if not sorted_bam:
            print('sambamba failed for '+star_bam)
            #log  
            error=sampleid+'\tSambamba failed'
            progresslog.write(error+'\n')
            continue

        #run stringtie assembly
        st=stieobj.perform_assembly(sorted_bam)
        #delete sorted bam
        if delete_bam:
            pe.delete_file(sorted_bam)
        if not st:
            print('Stringtie failed')
            #log
            error=sampleid+'\tStringTie failed'
            progresslog.write(error+'\n')
            continue
    if runquant:
        status=salmon.perform_quant(sobj)
        if not status:
            print('Salmon failed')
            #log
            error=sampleid+'\tSalmon failed'
            progresslog.write(error+'\n')
            continue
        #remove fastq
        sobj.delete_fastq()

    progresslog.write(sampleid+'\tcompleted'+'\n') 

#close log
progresslog.close()

