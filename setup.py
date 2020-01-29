#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 15:24:54 2019

@author: usingh
"""

import setuptools
import os
import sys



#exit if python 2
if sys.version_info.major != 3:
    raise EnvironmentError("""pyrpipe requires python 3.5 or higher.Please upgrade your python.""")
    
#read description
with open("README.md", "r") as fh:
    long_description = fh.read()

#read version info
cwd =os.path.abspath(os.path.dirname("__file__"))
version = {}
with open(os.path.join(cwd, "pyrpipe", "version.py")) as fp:
    exec(fp.read(), version)
version = version["__version__"]

if version is None:
    print("Error: version is missing. Exiting...", file=sys.stderr)
    sys.exit(1)



setuptools.setup(
    
    name="pyrpipe",
    version=version,
    author="Urminder Singh",
    author_email="usingh@iastate.edu",
    description="pyrpipe",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/urmi-21/pyrpipe",
    packages=setuptools.find_packages(),
    include_package_data=True,
    scripts=['scripts/pyrpipe_diagnostic.py'],
    install_requires=[line.rstrip() for line in open("requirements.txt", "rt")],
    tests_require=["pytest"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        
        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',


        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'
    ],
    python_requires='>=3.6',
    
)
