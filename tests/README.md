# Test instructions

It is recommended to create a conda environment for pyrpipe.

## Build conda environment

1. Run: `bash build_test_env.sh`
2. Activate the conda environment: `conda activate pyrpipe_test`
3. Please install cufflinks version 2.2.1 from [here](http://cole-trapnell-lab.github.io/cufflinks/install/)
4. Run tests from pyrpipe's root directory: `py.test tests/`