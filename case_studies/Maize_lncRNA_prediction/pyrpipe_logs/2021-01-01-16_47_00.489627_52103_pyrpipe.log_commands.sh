#!/bin/bash 
prefetch -O maize_out/SRR765545 SRR765545
fasterq-dump -O maize_out/SRR765545 -o SRR765545.fastq -e 6 -f maize_out/SRR765545/SRR765545.sra
prefetch -O maize_out/SRR765545 SRR765545
fasterq-dump -O maize_out/SRR765545 -o SRR765545.fastq -e 6 -f maize_out/SRR765545/SRR765545.sra
trim_galore --cores 6 --paired -o maize_out/SRR765545 maize_out/SRR765545/SRR765545_1.fastq maize_out/SRR765545/SRR765545_2.fastq
STAR --runMode genomeGenerate --genomeDir maize_out/starindex --genomeFastaFiles maize_out/Zea_mays.B73_RefGen_v4.dna.toplevel.fa --runThreadN 4
STAR --runMode genomeGenerate --genomeDir maize_out/starindex --genomeFastaFiles maize_out/Zea_mays.B73_RefGen_v4.dna.toplevel.fa --runThreadN 3
STAR --runMode genomeGenerate --genomeDir maize_out/starindex --genomeFastaFiles maize_out/Zea_mays.B73_RefGen_v4.dna.toplevel.fa --runThreadN 3
STAR --genomeChrBinNbits 10 --runMode genomeGenerate --genomeDir maize_out/starindex --genomeFastaFiles maize_out/Zea_mays.B73_RefGen_v4.dna.toplevel.fa --runThreadN 3