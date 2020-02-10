#!/bin/bash 
set -e

conda env create --force -f conda_env.yml
eval "$(conda shell.bash hook)"
conda activate pyrpipe_human
#update pip in case it fails
pip3 install --user --upgrade pip
pip install cython
pip install -Iv mikado==1.2.4
pip install git+https://github.com/urmi-21/pyrpipe.git

echo "Environment created. To activate, use: conda activate pyrpipe_human"