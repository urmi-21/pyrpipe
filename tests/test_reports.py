#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 13:23:55 2020

@author: usingh
"""

from pyrpipe import pyrpipe_engine as pe
from pyrpipe import pyrpipe_utils as pu

def test_report():
    cmd="pyrpipe_diagnostic.py report -o tests/testout/testreport tests/test_files/pyrpipe_logs/2020-01-22-18_14_47_pyrpipe.log"
    st=pe.execute_command(cmd.split(),verbose=True,quiet=False,logs=False,objectid="",command_name="")
    assert st==True, "report failed"


def test_benchmark():
    cmd="pyrpipe_diagnostic.py benchmark -t tests/testout/bmeport tests/test_files/pyrpipe_logs/2020-01-22-18_14_47_pyrpipe.log"
    st=pe.execute_command(cmd.split(),verbose=True,quiet=False,logs=False,objectid="",command_name="")
    assert st==True, "benchmark failed"
    
def test_shell():
    cmd="pyrpipe_diagnostic.py shell -o tests/testout/shellcmds tests/test_files/pyrpipe_logs/2020-01-22-18_14_47_pyrpipe.log"
    st=pe.execute_command(cmd.split(),verbose=True,quiet=False,logs=False,objectid="",command_name="")
    assert st==True, "shell failed"
    
    
def test_multiqc():
    cmd="pyrpipe_diagnostic.py multiqc -o tests/testout/mqcreport -t tests/testout/mqctmp tests/test_files/pyrpipe_logs/2020-01-22-18_14_47_pyrpipe.log"
    st=pe.execute_command(cmd.split(),verbose=True,quiet=False,logs=False,objectid="",command_name="")
    assert st==True, "shell failed"