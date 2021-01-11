import sys,os,json,yaml
from pyrpipe import pyrpipe_utils as pu
from pyrpipe import pyrpipe_engine as pe

"""
Download raw data and submit jobs
1: GTEx metadata file
2: out dir
3: slurm job name
Example:
python driver_GTEx.py adipose_metadata adipose_out adipose
"""

infile=sys.argv[1]
outdir=sys.argv[2]
jobid=sys.argv[3]
pu.mkdir(outdir)
cwd=os.getcwd()
#read config file
with open('config.yaml') as cf: config = yaml.load(cf, Loader=yaml.FullLoader)
batch_size=config['batch_size']
profile=config['gen3_profile']
threads=config['num_parallel']
job_time=str(config['slurm_job_time'])
dtn_ssh=config['dtn_ssh']
conda_profile=config['conda_profile']
conda_env=config['conda_env']
email=config['email']
pipeline=config['pipeline_type']

#load manifest file
print('reading manifest')
gtex_manifest = 'manifest.json'
with open(gtex_manifest, 'r') as fi:
    mdata = json.load(fi)
file_dict={}
for d in mdata:
    file_dict[d["file_name"]]=d


def extract_manifest(gtexids,outfile):
    result=[]
    for g in gtexids: result.append(file_dict[g+'.Aligned.sortedByCoord.out.patched.md.bam'])
    if not len(result) == len(gtexids):
        print('Some files not in manifest')
    #write manifest
    with open(outfile, 'w') as fo:
        json.dump(result , fo, indent=4)
    #write ids to file these will be used by process script
    fo=open(outfile.replace('.json','.ids'),'w')
    fo.write('\n'.join(gtexids)+'\n')
    fo.close()
     
def download_gtex_bams(manifest_file,outdir):
    #load list of bam files
    with open(manifest_file, 'r') as fi:
        thisdata = json.load(fi)
    flist=[]
    #check existing files
    for d in thisdata:
        f=d["file_name"] 
        gid=f.split('.Aligned')[0]
        outfile=os.path.join(outdir,gid,f)
        #if pu.check_files_exist(outfile) and pu.get_mdf(outfile)==d["md5sum"]:
        if pu.check_files_exist(outfile):
            print("Outfile {} exists. Skipping...".format(outfile))
            #copy it back to out dir
            os.rename(outfile,os.path.join(outdir,f))
        flist.append(d["file_name"]) 

    cmd='gen3-client download-multiple --profile={} --manifest={} --download-path={} --protocol=s3 --numparallel={} --skip-completed --no-prompt'.format(profile,m,outdir,threads)
    cdcmd='cd {}'.format(cwd)
    sshcmd=dtn_ssh+" '{}; {}'".format(cdcmd,cmd)
    out=pe.get_shell_output(sshcmd,verbose=True)
    
    #move the files
    for f in flist:
        source=os.path.join(outdir,f)
        gid=f.split('.Aligned')[0]
        destdir=os.path.join(outdir,gid)
        pu.mkdir(destdir)
        dest=os.path.join(destdir,f)
        #print('Moving {}-->{}'.format(source,dest))
        os.rename(source,dest)

#modify this as needed
def submit_slurm_job(outfile,jobname,time,command,email,conda_profile,conda_env):
    line='#!/bin/bash'
    line=line+'\n#SBATCH -J '+jobname
    line+='\n#SBATCH -N 1'
    line+='\n#SBATCH -p RM'
    line+='\n#SBATCH --ntasks-per-node 28'
    line+='\n#SBATCH -t '+time
    line+='\n#SBATCH -C EGRESS'
    line+='\n#SBATCH --mail-user='+email
    line+='\n#SBATCH --mail-type=ALL'
    line+='\n\n#load required modules'
    line+='\nsource '+conda_profile
    line+='\nsource ~/.bashrc'
    line+='\nconda activate '+conda_env

    line+='\n\n# run commands'
    cmds='\n'.join(command)
    line+='\n'+cmds+'\n'

    f=open(outfile,'w')
    f.write(line)
    f.close()

    #submiti
    print('Submitting job {}'.format(outfile))
    pe.execute_command(['sbatch',outfile])

###################################################################################################
#read input gtex ids
with open(infile) as fi:
    gtex_md=fi.read().splitlines()
gids=[]
for l in gtex_md:
    gids.append(l.split('\t')[0])

#split in batch size
f = lambda L, n=batch_size: [L[i:i+n] for i in range(0, len(L), n)]
print('write manifest')
#create manifest
manifest_list=[]
i=0
for l in f(gids):
    i+=1
    mout=os.path.join(outdir,'batch_'+str(i)+'.json')
    extract_manifest(l,mout)
    manifest_list.append(mout)
print('start download')
#download data in each manifest
for m in manifest_list:
    download_gtex_bams(m,outdir)
    #submit job to slurm
    outfile=m.replace('.json','.slurm')
    jobname=jobid+'_'+pu.get_filename(m).split('.json')[0]
    commands=[]
    commands.append('cd '+cwd)
    commands.append('python processGTEX.py '+m.replace('.json','.ids')+' '+pipeline)

    submit_slurm_job(outfile,jobname,job_time,commands,email,conda_profile,conda_env)


















