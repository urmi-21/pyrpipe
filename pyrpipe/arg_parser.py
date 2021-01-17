#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 18:17:58 2020

@author: usingh
"""

import argparse
import pyrpipe.version

ver=pyrpipe.version.__version__
parser = argparse.ArgumentParser(description='pyrpipe: A lightweight python package for RNA-Seq workflows (version {})'.format(ver),
usage="""
pyrpipe [<pyrpipe options>] --in <script.py> [<script options>]
OR
python <script.py> [<pyrpipe options>] [<script options>]

use pyrpipe_diagnostic command for reports, and tests and installation of RNA-Seq tools
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
parser.add_argument("--multiqc", help="Autorun MultiQC after execution\nDefault: False",default=False,dest='multiqc', action='store_true')

parser.add_argument("--version", help="Print version information and exit",default=False,dest='versioninfo', action='store_true')

#parser.add_argument('infile', help='The input python script',action="store",nargs="?")
required_input = parser.add_argument_group('Required input file if invoking via pyrpipe command')
required_input.add_argument('--in', help='The input python script', required=False,dest='infile')
