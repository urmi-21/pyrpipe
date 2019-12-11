#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 11:16:40 2019

@author: usingh
"""

import sys
import os
import argparse

"""
Valid options
help: print help

report: compile a full report of the analysis
shell: compile a shell script of all executed commands
failed: compile a list of failed commands
passed: compile a list of passed commands
benchmark: generate benchmarks from logs
"""

parser = argparse.ArgumentParser(prog='pyrpipe_diagnostic',add_help=True)
parser.add_argument("--verbose", help="increase output verbosity",action="store_true")
parser.add_argument("--report", help="increase output verbosity",action="store_true")
parser.add_argument("--shell", help="increase output verbosity",action="store_true")
parser.add_argument("--failed", help="increase output verbosity",action="store_true")
parser.add_argument("--passed", help="increase output verbosity",action="store_true")
parser.add_argument("--benchmark", help="increase output verbosity",action="store_true")
parser.add_argument("--logfile", help="increase output verbosity",action="store")







args = parser.parse_args()


###start processing
vFlag=args.verbose

if vFlag:
    print("pyrpipe diagnostic utility")


if args.report:
    print("Generating report")
