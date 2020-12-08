#!/bin/bash
## 2.1 indexing genome
bowtie2-build ref.fa example >index.stdout 2>index.stderr

## 2.2 alignment
for i in ./demultiplexed-barcode*.fastq.gz
 do
 echo $i;
 bowtie2 -p 8 --very-fast-local -x example -U $i -S ./${i##*/}.sam 2> ./${i##*/}.bowtie2.stdout;
done

## -p: number of processors to be used in alignment
## --very-fast-local: same as -D 5 -R 1 -N 0 -L 25 -i S,1,2.00
## -U: Comma-separated list of files containing unpaired reads to be aligned, e.g. lane1.fq,lane2.fq,lane3.fq,lane4.fq.
## -S: File to write SAM alignments to
