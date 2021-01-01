#!/bin/bash 
STAR --runMode genomeGenerate --genomeDir maize_out/starindex --genomeFastaFiles maize_out/Zea_mays.B73_RefGen_v4.dna.toplevel.fa --runThreadN 4
STAR --runMode genomeGenerate --genomeDir maize_out/starindex --genomeFastaFiles maize_out/Zea_mays.B73_RefGen_v4.dna.toplevel.fa --runThreadN 3
STAR --runMode genomeGenerate --genomeDir maize_out/starindex --genomeFastaFiles maize_out/Zea_mays.B73_RefGen_v4.dna.toplevel.fa --runThreadN 3
STAR --genomeChrBinNbits 10 --runMode genomeGenerate --genomeDir maize_out/starindex --genomeFastaFiles maize_out/Zea_mays.B73_RefGen_v4.dna.toplevel.fa --runThreadN 3