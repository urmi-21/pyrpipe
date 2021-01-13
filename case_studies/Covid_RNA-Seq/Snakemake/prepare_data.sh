#!/bin/bash
set -e

mkdir -p human_data
cd human_data

echo "[1] Downloading Transcriptome"
#download human transcriptome
wget -q ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.transcripts.fa.gz
#human genome
echo "[2] Downloading Genome"
wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

#download covid seq
echo "[3] SARS-CoV-2 Genome"
wget -q 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets' -O datasets
chmod +x datasets
./datasets download virus genome taxon 2697049 --refseq --filename SARS2-reference.zip
unzip SARS2-reference.zip
cp ncbi_dataset/data/genomic.fna sars.fna
rm -r SARS2-reference.zip ncbi_dataset 
#download viral decoys


echo "[4] Unzipping"
gunzip -f *.gz

echo "[5] Creating decoys"
#create decoy list
grep ">" GCA_000001405.15_GRCh38_no_alt_analysis_set.fna | cut -d " " -f 1 | tr -d ">" > decoys.txt
grep ">" sars.fna | cut -d " " -f 1 | tr -d ">" >> decoys.txt

#combine transcriptome and decoy fasta files
cat gencode.v36.transcripts.fa GCA_000001405.15_GRCh38_no_alt_analysis_set.fna sars.fna > human_tr_gen_decoy.fasta

echo "[6] Cleanup"
#cleanup
rm gencode.v36.transcripts.fa GCA_000001405.15_GRCh38_no_alt_analysis_set.fna 

#echo "[7] Creating Salmon index (this can take several minutes)"
#create salmon index
#time salmon index -t human_tr_gen_decoy.fasta -d decoys.txt -p 15 -i salmon_index

echo "[7] Done"
