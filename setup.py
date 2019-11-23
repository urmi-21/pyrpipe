#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 15:24:54 2019

@author: usingh
"""

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyrpipe", # Replace with your own username
    version="0.0.1",
    author="Urminder Singh",
    author_email="usingh@iastate.edu",
    description="pyrpipe",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/urmi-21/pyrpipe",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
    ],
    python_requires='>=3.6',
)
