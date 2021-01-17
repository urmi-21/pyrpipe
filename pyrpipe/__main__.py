#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 11:22:06 2020

@author: usingh

This is the main entry point into pyrpipe when invoked with pyrpipe comand
"""

import sys
import os
import pyrpipe.arg_parser
import pyrpipe.version


def main():
    ver=pyrpipe.version.__version__
    args, unknownargs = pyrpipe.arg_parser.parser.parse_known_args()   
    
    if args.versioninfo:
        print("pyrpipe version {}".format(ver))
        sys.exit(0)
           
    infile=args.infile    
    if not infile:
        pyrpipe.arg_parser.parser.print_help()
        sys.exit(1)
        
    args_list=args_to_list(args)
    
    cmd=['python',infile]+args_list+unknownargs
    
    #print('CMD',' '.join(cmd))
    
    os.system(' '.join(cmd))
    
    #perform auto reports    
    
    
def args_to_list(args):
    #threads
    threads=args.threads
    mem=args.mem
    paramdir=args.paramdir
    logsdir=args.logsdir
    dryrun=args.dryrun
    safemode=args.safemode
    nologs=args.nologs
    verbose=args.verbose
    force=args.force
    multiqc=args.multiqc
    
    args_list=[]
    
    if threads:
        args_list.append('--threads')
        args_list.append(threads)
    if mem:
        args_list.append('--max-memory')
        args_list.append(mem)
    if paramdir:
        args_list.append('--param-dir')
        args_list.append(paramdir)
    if logsdir:
        args_list.append('--logs-dir')
        args_list.append(logsdir)
    
    if verbose:
        args_list.append('--verbose')
    if dryrun:
        args_list.append('--dry-run')
    if force:
        args_list.append('--force')
    if safemode:
        args_list.append('--safe-mode')
    if nologs:
        args_list.append('--no-logs')
    if multiqc:
        args_list.append('--multiqc')
    
    #args_list.append('--build-tools')   
    
    return args_list
    
    
if __name__ == '__main__':
    main()


    
    
