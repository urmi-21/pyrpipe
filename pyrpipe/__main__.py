#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 11:22:06 2020

@author: usingh

This is the main entry point into pyrpipe
"""


import argparse
import sys
import os
import json
import pyrpipe.version




def main():
    ver=pyrpipe.version.__version__
    parser = argparse.ArgumentParser(description='pyrpipe: A lightweight python package for RNA-Seq workflows (version {})'.format(ver),
    usage="""
    pyrpipe [<options>] <infile.py>
    """)
    parser.add_argument("--threads", help="Num processes/threads to use\nDefault:mp.cpu_count()")
    parser.add_argument('--max-memory', help='Max memory to use (in GB)\ndefault: None',action="store",dest='mem')
    parser.add_argument("--verbose", help="Print pyrpipe_engine's stdout and stderr\nDefault: False",default=False,dest='verbose', action='store_true')
    parser.add_argument("--dry-run", help="Only print pyrpipe's commands and not execute anything through pyrpipe_engine module\nDefault: False",default=False,dest='dryrun', action='store_true')
    parser.add_argument("--force", help="Force execution of commands if their target files already exist\nDefault: False",default=False,dest='force', action='store_true')
    parser.add_argument("--safe-mode", help="Disable file deletions through pyrpipe_engine module\nDefault: False",default=False,dest='safemode', action='store_true')
    parser.add_argument("--no-logs", help="Disable pyrpipe logs\nDefault: False",default=False,dest='nologs', action='store_true')
    parser.add_argument("--param-dir", help="Directory containing parameter yaml files\nDefault: ./params",dest='paramdir',default='params')
    parser.add_argument("--logs-dir", help="Directory for saving log files\nDefault: ./pyrpipe_logs",dest='logsdir',default='pyrpipe_logs')
    
    parser.add_argument("--build-tools", help="Install the tools required by pyrpipe via  bioconda.",dest='build_conda',default=False,action='store_true')
    parser.add_argument("--test-sratools", help="Test if NCBI SRA-Tools is woking.",dest='test_sra',default=False,action='store_true')

    parser.add_argument("--version", help="Print version information and exit",default=False,dest='versioninfo', action='store_true')
    
    parser.add_argument('infile', help='The input python script',action="store",nargs="?")
    
    args = parser.parse_args()
    
    
    
    if args.versioninfo:
        print("pyrpipe version {}".format(ver))
        sys.exit(0)
           
    infile=args.infile    
    if not infile:
        parser.print_help()
        sys.exit(1)
    
    #threads
    procs=args.threads
    mem=args.mem
    dryrun=args.dryrun
    safemode=args.safemode
    paramdir=args.paramdir
    logsdir=args.logsdir
    nologs=args.nologs
    verbose=args.verbose
    force=args.force
    
    
        
    
    
    #call main program
    caller(procs,mem,dryrun,safemode,paramdir,logsdir,nologs,verbose,force,infile)
    
    
    
if __name__ == '__main__':
    main()

def caller(procs,mem,dryrun,safemode,paramdir,logsdir,nologs,verbose,force,infile):
    #write pyrpipe configuration
    #everything saved as str
    conf={}
    conf['threads']=procs
    conf['memory']=mem
    conf['params_dir']=paramdir
    conf['dry']=dryrun
    conf['force']=force
    conf['safe']=safemode
    conf['logs_dir']=logsdir
    conf['logging']= not nologs
    conf['verbose']=verbose
    
    with open('pyrpipe.conf', 'w') as outfile:
        json.dump(conf, outfile,indent=4)
    
    #execute
    cmd=['python',infile]
    os.system(' '.join(cmd))
    
    
    
