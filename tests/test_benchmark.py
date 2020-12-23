#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 22:27:28 2020

@author: usingh
"""
from pyrpipe import benchmark as bm

logfile='tests/test_files/pyrpipe_logs/2020-01-22-18_14_47_pyrpipe.log'
envlogfile='tests/test_files/pyrpipe_logs/2020-01-22-18_14_47_pyrpipeENV.log'
outdir='tests/testout/testreport'

def test1():
    ob=bm.Benchmark(logfile,envlogfile,out_dir=outdir)
    #generate benchmarks
    ob.plot_time_perobject()
    ob.plot_time_perprogram()