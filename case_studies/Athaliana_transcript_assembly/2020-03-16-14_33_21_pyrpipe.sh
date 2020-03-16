#!/bin/bash 
fasterq-dump -O athal_out/SRR976159 -o SRR976159.fastq -e 8 -f SRR976159
fasterq-dump -O athal_out/SRR978411 -o SRR978411.fastq -e 8 -f SRR978411
fasterq-dump -O athal_out/SRR971778 -o SRR971778.fastq -e 8 -f SRR971778
bbduk.sh in=athal_out/SRR976159/SRR976159_1.fastq in2=athal_out/SRR976159/SRR976159_2.fastq out=athal_out/SRR976159/SRR976159_1_bbduk.fastq out2=athal_out/SRR976159/SRR976159_2_bbduk.fastq threads=8 ktrim=r k=23 mink=11 qtrim='rl' trimq=10 ref=adapters2.fa -Xmx2g
bbduk.sh in=athal_out/SRR978411/SRR978411_1.fastq in2=athal_out/SRR978411/SRR978411_2.fastq out=athal_out/SRR978411/SRR978411_1_bbduk.fastq out2=athal_out/SRR978411/SRR978411_2_bbduk.fastq threads=8 ktrim=r k=23 mink=11 qtrim='rl' trimq=10 ref=adapters2.fa -Xmx2g
bbduk.sh in=athal_out/SRR971778/SRR971778_1.fastq in2=athal_out/SRR971778/SRR971778_2.fastq out=athal_out/SRR971778/SRR971778_1_bbduk.fastq out2=athal_out/SRR971778/SRR971778_2_bbduk.fastq threads=8 ktrim=r k=23 mink=11 qtrim='rl' trimq=10 ref=adapters2.fa -Xmx2g
hisat2 -1 athal_out/SRR976159/SRR976159_1_bbduk.fastq -2 athal_out/SRR976159/SRR976159_2_bbduk.fastq -S athal_out/SRR976159/SRR976159_hisat2.sam -p 10 -x athal_out/athalIndex/athalInd --dta-cufflinks
hisat2 -1 athal_out/SRR978411/SRR978411_1_bbduk.fastq -2 athal_out/SRR978411/SRR978411_2_bbduk.fastq -S athal_out/SRR978411/SRR978411_hisat2.sam -p 10 -x athal_out/athalIndex/athalInd --dta-cufflinks
hisat2 -1 athal_out/SRR971778/SRR971778_1_bbduk.fastq -2 athal_out/SRR971778/SRR971778_2_bbduk.fastq -S athal_out/SRR971778/SRR971778_hisat2.sam -p 10 -x athal_out/athalIndex/athalInd --dta-cufflinks
samtools view -o athal_out/SRR976159/SRR976159_hisat2.bam -@ 6 -b athal_out/SRR976159/SRR976159_hisat2.sam
samtools sort -o athal_out/SRR976159/SRR976159_hisat2_sorted.bam -@ 6 athal_out/SRR976159/SRR976159_hisat2.bam
samtools view -o athal_out/SRR978411/SRR978411_hisat2.bam -@ 6 -b athal_out/SRR978411/SRR978411_hisat2.sam
samtools sort -o athal_out/SRR978411/SRR978411_hisat2_sorted.bam -@ 6 athal_out/SRR978411/SRR978411_hisat2.bam
samtools view -o athal_out/SRR971778/SRR971778_hisat2.bam -@ 6 -b athal_out/SRR971778/SRR971778_hisat2.sam
samtools sort -o athal_out/SRR971778/SRR971778_hisat2_sorted.bam -@ 6 athal_out/SRR971778/SRR971778_hisat2.bam
stringtie -o athal_out/SRR976159/SRR976159_hisat2_sorted_stringtie.gtf -p 8 -G athal_out/Arabidopsis_thaliana.TAIR10.45.gtf athal_out/SRR976159/SRR976159_hisat2_sorted.bam
stringtie -o athal_out/SRR978411/SRR978411_hisat2_sorted_stringtie.gtf -p 8 -G athal_out/Arabidopsis_thaliana.TAIR10.45.gtf athal_out/SRR978411/SRR978411_hisat2_sorted.bam
stringtie -o athal_out/SRR971778/SRR971778_hisat2_sorted_stringtie.gtf -p 8 -G athal_out/Arabidopsis_thaliana.TAIR10.45.gtf athal_out/SRR971778/SRR971778_hisat2_sorted.bam