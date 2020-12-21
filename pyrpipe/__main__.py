#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 11:22:06 2020

@author: usingh

This is the main entry point into pyrpipe
"""

import sys
import os
import json
import pyrpipe.arg_parser
import pyrpipe.version




def main():
    ver=pyrpipe.version.__version__
    #args = parser.parse_args()
    args, unknownargs = pyrpipe.arg_parser.parser.parse_known_args()
    
    
    
    if args.versioninfo:
        print("pyrpipe version {}".format(ver))
        sys.exit(0)
           
    infile=args.infile    
    if not infile:
        pyrpipe.arg_parser.parser.print_help()
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
    caller(procs,mem,dryrun,safemode,paramdir,logsdir,nologs,verbose,force,infile,unknownargs)
    
    
    
if __name__ == '__main__':
    main()

def caller(procs,mem,dryrun,safemode,paramdir,logsdir,nologs,verbose,force,infile,unknownargs):
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
    cmd=['python',infile]+unknownargs
    os.system(' '.join(cmd))
    
    
    
