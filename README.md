# snpGBS


Previously known as `HOMEBREW`, developed by **Rudiger Brauning**.


## Overview

Existing SNP calling pipelines often come with built-in filtering processes, which can introduce **systematic bias**. Each method has its own algorithm to identify and define SNPs, it can lead to **moderate to large variations in final outputs**.
A **simple**, **intuitive** and **consistent** bioinformatics workflow is thus needed for developing new analytical methods.

Here we present `snpGBS`, a simple **three-step approach to identify SNPs from GBS data**:



### Step One: Demultiplexing

We use `cutadapt` to demultiplex the raw GBS data (i.e. **.fastq** or **.fastq.gz** file). More information about `cutadapt`: <https://cutadapt.readthedocs.io/en/stable/guide.html#demultiplexing>

```bash
## 1.1 trimming common adapter
cutadapt -j 8 -a common_adapter=AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG -o example.trimmed.fastq example.fastq >01.trimmed.stdout 2>01.trimmed.stderr

## 1.2 demultiplexing
cutadapt -j 8 -e 0 --no-indels -g file:barcodes.fasta -o "demultiplexed_{name}.fastq.gz" example.trimmed.fastq >01.demultiplexed.stdout 2>01.demultiplexed.stderr
```



### Step Two: Mapping

We use `bowtie2` to align and map GBS reads. More information about `bowtie2`: <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>

```bash
## 2.1 indexing genome
bowtie2-build ref.fa example >index.stdout 2>index.stderr

## 2.2 alignment
for i in ./demultiplexed_barcode*.fastq.gz
 do
    echo $i;
    bowtie2 --very-fast-local -x example -U $i -S ./${i##*/}.sam 2>./${i##*/}.bowtie2.stdout;
done
```


### Step Three: SNP Calling

We use `bcftools-mpileup` to identify SNPs. More information about `bcftools-mpileup`: <http://www.htslib.org/doc/bcftools.html#mpileup>

```bash
## convert SAM to BAM
## 3.1 convert SAM to BAM
for i in *.sam;
  do
    echo $i;
    samtools view -bS $i > "${i%.sam}.bam";
done

## 3.2 sort bam
for i in *.bam;
  do
    echo $i;
    samtools sort $i -o "${i%.bam}.sorted.bam";
done

## 3.3 create bamlist
for i in *.sorted.bam;
  do
    echo $i;
done > bamlist;


## 3.4 calling SNPs
bcftools mpileup -I -Ou -f ref.fa -b bamlist -a AD | bcftools call -cv - | bcftools view -M2 - >example.vcf
```


### Required Input

#### Demultiplexing

- **Raw GBS data** (e.g. `example.fastq`)

- **Barcode sequences** (e.g. `barcodes.fasta`)

  `createBarcodeFASTA.sh` is provided to convert `barcodes.txt` to `barcodes.fasta` as below

  ```bash
  #!/bin/bash
  filename='barcodes.txt'
  n=0
  while read line; do
   echo ">barcode$n"
   echo "^$line"
   n=$((n+1))
  done < $filename
  ```

#### Mapping

- **Reference genome** (e.g.`ref.fa`)

#### SNP Calling

- **Nil**


### Example

We help users with testing `snpGBS`, we put together an example with the following files

- **Raw GBS data**: `example.fastq`can be found in <https://figshare.com/articles/dataset/snpGBS/13591274>

- **Barcode sequences**: `barcodes.txt` and `barcodes.fasta` are stored in <https://github.com/AgResearch/snpGBS/tree/main/example/datasets>

- **Reference Genome**: `ref.fa` can also be found in <https://figshare.com/articles/dataset/snpGBS/13591274>

### Output

Here's a list of expected outputs

- **Demultiplexed FASTQ files**:  `demultiplexed-barcode*.fastq.gz` in https://github.com/AgResearch/snpGBS/tree/main/example/demultiplexed_fastq

- **Mapped, sorted and indexed BAM files**: `demultiplexed-barcode*.fastq.gz.sorted.bam(.bai)` in https://github.com/AgResearch/snpGBS/tree/main/example/mapping

- **VCF file**: `example.vcf` stores all the SNP information, it can be downloaded from <https://figshare.com/articles/dataset/snpGBS/13591274>
