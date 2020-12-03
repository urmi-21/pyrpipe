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
import subprocess




def main():
    ver=pyrpipe.version.__version__
    print('pyrpipe version',ver)
    #set pyrpipe configuration
    conf={}
    if '--dry-run' in sys.argv:
        conf['dry']=True
    else:
        conf['dry']=False
        
    #write file
    with open('pyrpipe.conf', 'w') as outfile:
        json.dump(conf, outfile)
    
    #execute
    print('executing input file',sys.argv[1])
    cmd=['python',sys.argv[1]]
    print('CMD',cmd)
    #result = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
    #stdout,stderr = result.communicate()
    
    #print(result.returncode,stdout,stderr)
    os.system(' '.join(cmd))
    
    print('Ending CMD',cmd)
    
    
    
    
if __name__ == '__main__':
    main()