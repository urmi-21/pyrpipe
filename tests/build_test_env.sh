#!/bin/bash 
set -e

conda env create --force -f test_environment.yml
eval "$(conda shell.bash hook)"
conda activate pyrpipe_test
#update pip in case it fails
pip install --user --upgrade pip
pip install cython
pip install -Iv mikado==1.2.4

echo "Environment created. To activate, use: conda activate pyrpipe_test"