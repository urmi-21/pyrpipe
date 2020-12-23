#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 17:57:29 2020

@author: usingh
"""

from pyrpipe import pyrpipe_utils as pu

def test_time():
    print(pu.get_timestamp(shorten=True))
    print(pu.get_timestamp(shorten=False))
    
def test_file_path():
    assert pu.check_files_exist('tests/test_utils.py')==True, "Check file failed"
    assert pu.check_files_exist('tests/test_utils.py','tests/test_qc.py')==True, "Check file failed"
    assert pu.check_files_exist('tests/XXX.py')==False, "Check file failed"
    
def test_path():
    assert pu.check_paths_exist('tests/')==True, "Check path failed"
    assert pu.check_paths_exist('tests/XXXXtest_utils.py','tests/test_qc.py')==False, "Check path failed"
    
def test_bytetoread():
    assert pu.byte_to_readable(12312)=='12.0 KB', 'byte to readable failed'
    
def test_java_args():
    res=['A=3', 'B=22', '-C']
    assert pu.parse_java_args(['A','B','-C'], {"A": "3", "B": "22","-C":""})==res,'java args failed'
    
def test_linux_args():
    res=['-O', './test', 'IN1', 'IN2']
    assert pu.parse_unix_args(['-O','-t','-q'], {"-O": "./test", "Attr2": "XX","--":("IN1","IN2")})==res, 'linux args failed'
    
def test_get_dir():
    assert pu.get_file_directory('aaa/bb')=='aaa', 'get_dir failed'
    assert pu.get_file_directory('bbb/aaa/bb')=='bbb/aaa', 'get_dir failed'
    
def test_filename():
    assert pu.get_filename('aaa/bbb.tss')=='bbb.tss','get filename failed'
    
def test_ext():
    assert pu.get_fileext('aaa/bbb.tss')=='.tss','get fileext failed'
    
def test_mkdir():
    assert pu.mkdir('tests/testout/testdira')==True,'mkdir failed'
    assert pu.mkdir('tests/testout/testdira')==False,'mkdir failed'
    
def test_union():
    pu.get_union([1,2,3,4],[1,2,3,5]).sort()==[1,2,3,4,5].sort(),'union failed'
    
