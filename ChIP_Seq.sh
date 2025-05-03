#!/bin/bash

# Install SRA toolkit
http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz
tar -vxzf /home/nuria/Documents/PGB/Chip/Olig2/sratoolkit.current-mac64.tar.gz
export PATH=$PATH:/home/nuria/Documents/PGB/Chip/Olig2/sratoolkit.3.0.7-mac64/bin

# Download input and peak raw data
fastq-dump -v --split-files --gzip --outdir SRR1583894
fastq-dump -v --split-files --gzip --outdir SRR1583890

# Quality control analysis using FastQC
fastqc SRR1583894_1.fastq.gz
fastqc SRR1583890_1.fastq.gz

# Summary of the samples quality control using MultiQC
multiqc . *

#Trim low-quality 3' ends from reads (Phred score < 30) and discard reads
shorter than 20bp.
cutadapt -m 20 -q 30 -o SRR1583890_1_trim.fastq.gz SRR1583890_1.fastq.gz

# New quality control to analyze improvements
fastqc SRR1583890_1_trim.fastq.gz
multiqc . *

# Map sample reads with the reference genome using bwa
bwa mem /home/nuria/Documents/PGB/Chip/mousegenome SRR1583894_1_trim.fastq >
SRR1583894.sam
bwa mem /home/nuria/Documents/PGB/Chip/mousegenome SRR1583890_1_trim.fastq >
SRR1583890.sam

# Sort and compress .sam files using samtools
samtools sort -O bam SRR1583894.sam > SRR1583894.bam
samtools sort -O bam SRR1583890.sam > SRR1583890.bam

# Indexing .bam files
samtools index SRR1583894.bam
samtools index SRR1583890.bam

# Peak calling using MACS2
macs2 callpeak -g mm -f BAM -t SRR1583890.bam -c SRR1583894.bam --bw 200 -q 0.1
-outdir . -n Olig2

# How many peaks do we have?
wc -l Olig2_peaks.narrowPeak
# Output: 24256 peaks

# Resize peaks to equal length from midpoint
awk '{print $1,$2+int(($3-$2)/2)-250,$2+int(($3-$2)/2)+250}' Olig2_peaks.narrowPeak | tr " " "\t" | sort -k1,1 -k2,2n > Peaks.resized.bed

# Obtain sequences from the reference genome using bedtools
fastaFromBed -fi /home/nuria/Documents/PGB/Chip/mousegenome -bed Peaks.resized.bed > Peaks.resized.bed.fa

# Install MEME-CHIP
tar zxf meme-5.5.4.tar.gz
cd meme-5.5.4
./configure --prefix=$HOME/meme --enable-build-libxml2 --enable-build-libxslt
make
make test
make install

# MEME-CHIP
meme-chip -oc Olig2_meme --db jaspar.meme Peaks.resized.bed.fa
