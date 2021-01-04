#!/bin/bash

mkdir -p refdata
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-49/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz -O refdata/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-49/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.49.gff3.gz -O refdata/Saccharomyces_cerevisiae.R64-1-1.49.gff3.gz

cd refdata
gunzip -f *.gz


#run star index
mkdir -p yeast_index
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir yeast_index --genomeFastaFiles Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa --genomeSAindexNbases 10

