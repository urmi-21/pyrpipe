# Introduction
## Prediction of long non-coding RNAs (lncRNAs) in *Zea mays* by supplementing pyrpipe with a third-party tool

We downloaded RNA-Seq data, quality filtered using Trim Galore, aligned reads using STAR and transcripts were assembled using StringTie.
Then, we used a third party tool, (PLncPRO), to predict lncRNAs in the assembled transcripts.

## To run the example:

1. Build the conda environment:

    `conda env create -f environment.yml`
    
2. Switch to the environment:

    `conda activate pyrpipe_examples`

3. Execute the code in `.ipynb` or `.py` file
