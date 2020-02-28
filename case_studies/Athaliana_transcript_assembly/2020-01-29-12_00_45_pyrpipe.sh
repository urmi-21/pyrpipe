#!/bin/bash 
prefetch -O athal_out/SRR976159 SRR976159
prefetch -O athal_out/SRR978411 SRR978411
prefetch -O athal_out/SRR971778 SRR971778
fasterq-dump -e 8 -f -t athal_out -O athal_out/SRR976159 -o SRR976159.fastq athal_out/SRR976159/SRR976159.sra
fasterq-dump -e 8 -f -t athal_out -O athal_out/SRR978411 -o SRR978411.fastq athal_out/SRR978411/SRR978411.sra
fasterq-dump -e 8 -f -t athal_out -O athal_out/SRR971778 -o SRR971778.fastq athal_out/SRR971778/SRR971778.sra
bbduk.sh ktrim=r k=23 mink=11 qtrim='rl' trimq=10 ref=adapters2.fa in=athal_out/SRR976159/SRR976159_1.fastq in2=athal_out/SRR976159/SRR976159_2.fastq out=athal_out/SRR976159/SRR976159_1_bbduk.fastq out2=athal_out/SRR976159/SRR976159_2_bbduk.fastq -Xmx2g
bbduk.sh ktrim=r k=23 mink=11 qtrim='rl' trimq=10 ref=adapters2.fa in=athal_out/SRR978411/SRR978411_1.fastq in2=athal_out/SRR978411/SRR978411_2.fastq out=athal_out/SRR978411/SRR978411_1_bbduk.fastq out2=athal_out/SRR978411/SRR978411_2_bbduk.fastq -Xmx2g
bbduk.sh ktrim=r k=23 mink=11 qtrim='rl' trimq=10 ref=adapters2.fa in=athal_out/SRR971778/SRR971778_1.fastq in2=athal_out/SRR971778/SRR971778_2.fastq out=athal_out/SRR971778/SRR971778_1_bbduk.fastq out2=athal_out/SRR971778/SRR971778_2_bbduk.fastq -Xmx2g
hisat2-build -p 8 -a -q athal_out/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa athal_out/athalIndex/athalInd
hisat2 --dta-cufflinks -p 10 -x athal_out/athalIndex/athalInd -1 athal_out/SRR976159/SRR976159_1_bbduk.fastq -2 athal_out/SRR976159/SRR976159_2_bbduk.fastq -S athal_out/SRR976159/SRR976159_hisat2.sam
hisat2 --dta-cufflinks -p 10 -x athal_out/athalIndex/athalInd -1 athal_out/SRR978411/SRR978411_1_bbduk.fastq -2 athal_out/SRR978411/SRR978411_2_bbduk.fastq -S athal_out/SRR978411/SRR978411_hisat2.sam
hisat2 --dta-cufflinks -p 10 -x athal_out/athalIndex/athalInd -1 athal_out/SRR971778/SRR971778_1_bbduk.fastq -2 athal_out/SRR971778/SRR971778_2_bbduk.fastq -S athal_out/SRR971778/SRR971778_hisat2.sam
samtools view -@ 8 -o athal_out/SRR976159/SRR976159_hisat2.bam -b athal_out/SRR976159/SRR976159_hisat2.sam
samtools sort -@ 8 -o athal_out/SRR976159/SRR976159_hisat2_sorted.bam athal_out/SRR976159/SRR976159_hisat2.bam
samtools view -@ 8 -o athal_out/SRR978411/SRR978411_hisat2.bam -b athal_out/SRR978411/SRR978411_hisat2.sam
samtools sort -@ 8 -o athal_out/SRR978411/SRR978411_hisat2_sorted.bam athal_out/SRR978411/SRR978411_hisat2.bam
samtools view -@ 8 -o athal_out/SRR971778/SRR971778_hisat2.bam -b athal_out/SRR971778/SRR971778_hisat2.sam
samtools sort -@ 8 -o athal_out/SRR971778/SRR971778_hisat2_sorted.bam athal_out/SRR971778/SRR971778_hisat2.bam
stringtie -G athal_out/Arabidopsis_thaliana.TAIR10.45.gtf -o athal_out/SRR976159/SRR976159_hisat2_sorted_stringtie.gtf athal_out/SRR976159/SRR976159_hisat2_sorted.bam
stringtie -G athal_out/Arabidopsis_thaliana.TAIR10.45.gtf -o athal_out/SRR978411/SRR978411_hisat2_sorted_stringtie.gtf athal_out/SRR978411/SRR978411_hisat2_sorted.bam
stringtie -G athal_out/Arabidopsis_thaliana.TAIR10.45.gtf -o athal_out/SRR971778/SRR971778_hisat2_sorted_stringtie.gtf athal_out/SRR971778/SRR971778_hisat2_sorted.bam