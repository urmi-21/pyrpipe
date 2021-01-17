#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 12:04:28 2019

@author: usingh
"""

import os
import datetime as dt
import sys
from colorama import Fore,Back,Style
import re
import hashlib


def print_boldred(text):
    """Print in bold red font
    
    Parameters
    ----------
    
    text: str
        text to print
        
    Returns:   None
    """
    #print (Colors.FAIL + Colors.BOLD+ text + Colors.ENDC)
    print_error(text)

def print_green(text):
    """Print in green font
    
    Parameters
    ----------
    
    text: str
        text to print
        
    Returns:   None
    """
    #print (Colors.OKGREEN + text + Colors.ENDC)
    print_success(text)

def print_blue(text):
    """Print in blue font
    
    Parameters
    ----------
    
    text: str
        text to print
        
    Returns:   None
    """
    #print (Colors.OKBLUE + text + Colors.ENDC) 
    print_message(text)

def print_yellow(text):
    """Print in yellow font
    
    Parameters
    ----------
    
    text: str
        text to print
        
    Returns:   None
    """
    #print (Colors.LightYellow + text + Colors.ENDC) 
    print_notification(text)
    

def pyrpipe_print(color,*args,stderr=False,**kwargs):
    if stderr:
        print(color,file=sys.stderr,end="",**kwargs)
        print(*args,file=sys.stderr,end="",**kwargs)
        print(Style.RESET_ALL,file=sys.stderr,**kwargs)
    else:
        print(color,file=sys.stdout,end="",**kwargs)
        print(*args,file=sys.stdout,end="",**kwargs)
        print(Style.RESET_ALL,file=sys.stdout,**kwargs)
    
def print_notification(*args):
    """Print message to stderr
    """
    #args=('[pyrpipe notification]',)+args
    pyrpipe_print(Fore.LIGHTYELLOW_EX,*args,stderr=True)

def print_message(*args):
    """Print message to stderr
    """
    pyrpipe_print(Fore.LIGHTCYAN_EX,*args,stderr=True)
    
def print_warning(*args):
    """Print message to stderr
    """
    args=('WARNING:',)+args
    pyrpipe_print(Fore.LIGHTMAGENTA_EX,*args,stderr=True)
    
def print_error(*args):
    """Print message to stderr
    """
    pyrpipe_print(Fore.LIGHTRED_EX,*args,stderr=True)
    
def print_success(*args):
    """Print message to stderr
    """
    pyrpipe_print(Fore.LIGHTGREEN_EX,*args,stderr=True)

######End color functions###################

def get_timestamp(shorten=False):
    """Function to return current timestamp.
    
    Parameters
    ----------
    
    shorten: bool
        return short version without space, dash and colons

    :return: timestamp as string
    :rtype: string
    """
    
    timestamp=str(dt.datetime.now()).split(".")[0].replace(" ","-")
    if shorten:
        timestamp=timestamp.replace("-","").replace(" ","").replace(":","")
    return timestamp
    

def check_paths_exist(*args):
    """Function to check if a directory exists.
    
    Parametrs
    ---------
    
    args: tuple
        a list of paths to check

    :return: return true only if all paths exist
    :rtype: bool
    """
    fail_flag=False
    for path in args:
        if not (os.path.exists(path) and os.path.isdir(path)):
            fail_flag=True
    if fail_flag==True:
        return False
    return True


def check_files_exist(*args):
    """Function to check if files exist.
    
    Parametrs
    ---------
    
    args: tuple
        a list of paths to check

    :return: return true only if all files exist
    :rtype: bool
    """
    for path in args:
        if not path:
            return False
        elif not os.path.isfile(path):
            return False
            
    return True


def check_hisatindex(index):
    """Function to check if hisat2 index is valid and exists.
    
    Parameters
    ----------
    
    index: str
        Path to the index 

    :return: Return true if index is valid
    :rtype: bool
    """
    if not index:
        return False
    if not check_files_exist(index+".1.ht2"):
        return check_files_exist(index+".1.ht2l")
    return True

def check_kallistoindex(index):
    """Function to check if kallisto index is valid and exists.
    
    Parameters
    ----------
    
    index: str
        Path to the index 

    :return: Return true if index is valid
    :rtype: bool
    """
    if not index:
        return False
    if not check_files_exist(index):
        return False
    return True

def check_salmonindex(index):
    """Function to check if salmon index is valid and exists.
    
    Parameters
    ----------
    
    index: str
        Path to the index 

    :return: Return true if index is valid
    :rtype: bool
    """
    if not index:
        return False
    if not check_paths_exist(index):
        return False
    return True


def check_starindex(index):
    """Function to check if star index is valid and exists.
    
    Parameters
    ----------
    
    index: str
        Path to the index 

    :return: Return true if index is valid
    :rtype: bool
    """
    if not index:
        return False
    
    if check_paths_exist(index):
        files_to_check=['chrLength.txt',
                      'chrNameLength.txt',
                      'chrName.txt',
                      'chrStart.txt',
                      'genomeParameters.txt',
                      'Genome']
        for f in files_to_check:
            if not check_files_exist(os.path.join(index,f)):
                return False
        return True
    
    return False


def check_bowtie2index(index):
    """Function to check if bowtie2 index is valid and exists.
    
    Parameters
    ----------
    
    index: str
        Path to the index 

    :return: Return true if index is valid
    :rtype: bool
    """
    if not index:
        return False
    if not check_files_exist(index+".1.bt2"):
        return check_files_exist(index+".1.bt2l")
    return True
    

def byte_to_readable(size_bytes):
    """Function to convert bytes to human readable format (MB,GB ...)
    
    Parameters
    ----------
    
    size_bytes: float
        size in bytes

    :return: Return size in human readable format
    :rtype: str
    """
    
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return "%3.1f %s" % (size_bytes, x)
        size_bytes /= 1024.0


def get_file_size(file_path):
    """Returns file size in human readable format
    
    Parameters
    ----------
    
    file_path: str
        Path to file

    :return: Return size in human readable format
    :rtype: string
    """
    
    if (check_files_exist(file_path)):
        file_info = os.stat(file_path)
        return byte_to_readable(file_info.st_size)
    

#TODO: override in case of empty list
def parse_java_args(valid_args_list,passed_args):
    """
    Function creates arguments to pass to java programs
    
    Parameters
    ----------
    
    valid_args_list: list
        list of valid arguments. Invalid arguments will be ignored
    passed_args: *dict
        keyword value argument list to be parsed

    :return: a list with command line arguments to be used with subprocess.popen
    :rtype: list
    
    Examples
    --------
    >>> parse_java_args(['A','B','-C'], {"A": "3", "B": "22","-C":""})
        ['A=3', 'B=22', '-C']
    """
    popen_args=[]
    special_args=["--"]
    positional_args=[]
    
    #empty list supplied consider all armunets valid
    if not valid_args_list:
        valid_args_list=list(passed_args.keys())
        #above command will also add specialArgs, remove those
        for x in special_args:
            if x in valid_args_list:
                valid_args_list.remove(x)
    
    for key, value in passed_args.items():
        #check if key is a valid argument
        if key in valid_args_list:
            #do not add emty parameters e.g. -q or -v
            if len(value)>0:
                popen_args.append(key+"="+value)
            else:
                popen_args.append(key)
        elif key in special_args:
            positional_args.extend(value)
        
        else:
            print_warning("Ignoring invalid argument '{0}: {1}'...".format(key, value))
    popen_args.extend(positional_args)
    return popen_args
    
    
def parse_unix_args(valid_args_list,passed_args):
    """Function creates command line arguments to pass to unix programs
    
    Parameters
    ----------
    
    valid_args_list: list
        list of valid arguments. Invalid arguments will be ignored
    passed_args: *dict
        keyword value argument list to be parsed
        
    :return: a list with command line arguments to be used with subprocess.popen
    :rtype: list
    
    Examples
    --------
    >>> parse_unix_args(['-O','-t','-q'], {"-O": "./test", "Attr2": "XX","--":("IN1","IN2")})
        Unknown argument Attr2 XX. ignoring...
        ['-O', './test', 'IN1', 'IN2']
    """
    popen_args=[]
    """
    Define some special arguments.
    -- to pass input which don't follow any flag e.g. mypro.sh input1 input2
    """
    special_args=["--"]
    positional_args=[]
    
    #empty list or None supplied consider all armunets valid
    if not valid_args_list:
        valid_args_list=list(passed_args.keys())
        #above command will also add specialArgs, remove those
        for x in special_args:
            if x in valid_args_list:
                valid_args_list.remove(x)
        
        
    for key, value in passed_args.items():
        #check if key is a valid argument
        if key in valid_args_list:
            popen_args.append(key)
            #do not add emty parameters e.g. -q or -v
            if len(value)>0:
                #if there are multiple values for a flag separated by space
                popen_args.extend(value.split(" "))
                        
                #popen_args.append(value)
        elif key in special_args:
            positional_args.extend(value)
        else:
            print_warning("Ignoring invalid argument '{0}: {1}'...".format(key, value))
    popen_args.extend(positional_args)
    return popen_args
    




def get_file_directory(file_path):
    """Returns directory of a file
    
    Parameters
    ----------
    
    file_path: str
        Path to file

    :return: directory os file_path
    :rtype: string
    """
    return os.path.split(file_path)[0]

def get_filename(file_path):
    """Returns filename with extension
    
    Parameters
    ----------
    
    file_path: str
        Path to file

    :return: filename
    :rtyep: string
    """
    return os.path.split(file_path)[1]

def get_fileext(file_path):
    """Returns file extension
    
    Parameters
    ----------
    
    file_path: str
        Path to file

    :return: file extension
    :rtype: string
    """
    return os.path.splitext(file_path)[1]

def get_file_basename(file_path):
    """Returns file basename without extension
    
    Parameters
    ----------
    
    file_path: str
        Path to file

    :return: file basename without extension
    :rtype: string
    """
    return os.path.splitext(get_filename(file_path))[0]
    


def mkdir(dir_path):
    """Create a directory

    :return: true is directory created
    :rtype: bool
    """
    try:
        os.makedirs(dir_path)
    except OSError:
        return False
    return True


def get_union(*args):
    """Return unioin of multiple input lists.
    """
    return list(set().union(*args)) 
    
def find_files(search_path,search_pattern,recursive=False,verbose=False):
    """
    Example: find_files(path,".fastq$")

    Parameters
    ----------
    search_path : TYPE
        DESCRIPTION.
    search_pattern : TYPE
        DESCRIPTION.
    recursive : TYPE, optional
        DESCRIPTION. The default is False.
    verbose : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    result : TYPE
        DESCRIPTION.
    
    Example
    --------
    find_files('path/to/directory','.*\.txt$',recursive=False,verbose=False)

    """
    
    result=[]
    if not search_path:
        search_path='./'
    #if search path is not valid
    if not check_paths_exist(search_path):
        return result
    
    pattern = re.compile(search_pattern)
    
    if recursive:
        for root, dirs, files in os.walk(search_path):
            for file in files:
                if bool(pattern.search(file)):
                    result.append(os.path.join(root,file))
        return result
    
    for file in os.listdir(search_path):
        if bool(pattern.search(file)):
            result.append(os.path.join(search_path,file))
    
    
    return result
             

def get_mdf(filename):
    """
    Compute and return md5checksum 

    Parameters
    ----------
    filename : str
        Path to input file

    Returns
    -------
    md5 : str
        MD5 checksum value

    """
    if not check_files_exist(filename):
        #print_boldred('Error calculating MD5. File {} not found'.format(filename))
        return None
    with open(filename,"rb") as f:
    # read contents
        data = f.read()    
        # pipe contents of the file through
    md5 = hashlib.md5(data).hexdigest()
    return md5
    


