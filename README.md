# snpGBS

`snpGBS` is a bioinformatics workflow to identify **single nucleotide polymorphism (SNP)** from **Genotyping-by-Sequencing (GBS)** data.

Previously known as `HOMEBREW`, developed by **Rudiger Brauning**.


## Overview

Existing SNP calling pipelines/software often come with built-in filtering processes, which can introduce **systematic bias**. Each different method has its own algorithm to identify and define SNPs, and it can lead to **moderate to large variations in final outputs**. A **simple**, **intuitive** and **consistent** bioinformatics workflow is thus needed for developing new analytical methods.

Here we present `snpGBS`, a **three-step approach to identify SNPs from GBS data**:



### Step One: Demultiplexing

We use `cutadapt` to demultiplex the raw GBS data (i.e. **.fastq** or **.fastq.gz** file)

More information about `cutadapt`: <https://cutadapt.readthedocs.io/en/stable/guide.html#demultiplexing>

```bash
cutadapt -j 8 -e 0 --no-indels -a common_adapter=AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG -g file:barcodes.fasta -o "demultiplexed_{name}.fastq.gz" example.fastq >01.demultiplexed.stdout 2>01.demultiplexed.stderr
```



### Step Two: Mapping

We use `bowtie2` to align and map GBS reads.

More information about `bowtie2`: <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>

```bash
## 2.1 indexing genome
bowtie2-build ref.fa example >index.stdout 2>index.stderr

## 2.2 alignment
for i in ./demultiplexed-barcode*.fastq.gz
 do
    echo $i;
    bowtie2 --very-fast-local -x example -U $i -S ./${i##*/}.sam 2>./${i##*/}.bowtie2.stdout;
done
```


### Step Three: SNP Calling

We use `bcftools-mpileup` to identify SNPs.

More information about `bcftools-mpileup`: <http://www.htslib.org/doc/bcftools.html#mpileup>

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

**To test `snpGBS`, we provide the following files in the `example` folder**

- `barcodes.txt` and `barcodes.fasta`


**`example.fastq`, `ref.fa` files (and the expected outputs) can be found in `/dataset/GBS_Sim/archive/snpGBS`.**


### Conda Environment

We built a `conda` **enviroment** located in `/home/kangj/conda-envs/homebrew`. You can activate it by

```bash
conda activate /home/kangj/conda-envs/homebrew/
```

It contains the following **packages**:

| Name                   | Version           | Build            | Channel    |
| ---------------------  | -----------       |-----------       | -----------|
| _libgcc_mutex          |   0.1             |    conda_forge   | conda-forge |
| _openmp_mutex          |   4.5             |          1_gnu   | conda-forge |
| bcftools               |   1.9             |     h68d8f2e_9   | bioconda    |
| bowtie2                |   2.4.2           | py38h1c8e9b9_1   | bioconda    |
| bzip2                  |   1.0.8           |     h7f98852_4   | conda-forge |
| c-ares                 |   1.17.1          |     h36c2ea0_0   | conda-forge |
| ca-certificates        |   2020.11.8       |     ha878542_0   | conda-forge |
| certifi                |   2020.11.8       | py38h578d9bd_0   | conda-forge |
| curl                   |   7.71.1          |     he644dc0_8   | conda-forge |
| cutadapt               |   3.0             | py38h0213d0e_0   | bioconda    |
| dnaio                  |   0.4.4           | py38h0213d0e_0   | bioconda    |
| gsl                    |   2.5             |     h294904e_1   | conda-forge |
| htslib                 |   1.9             |     h4da6232_3   | bioconda    |
| isa-l                  |   2.30.0          |     h36c2ea0_0   | conda-forge |
| krb5                   |   1.17.2          |     h926e7f8_0   | conda-forge |
| ld_impl_linux-64       |   2.35.1          |     hed1e6ac_0   | conda-forge |
| libblas                |   3.9.0           |     3_openblas   | conda-forge |
| libcblas               |   3.9.0           |     3_openblas   | conda-forge |
| libcurl                |   7.71.1          |     hcdd3856_8   | conda-forge |
| libdeflate             |   1.7             |     h36c2ea0_0   | conda-forge |
| libedit                |   3.1.20191231    |     he28a2e2_2   | conda-forge |
| libev                  |   4.33            |     h516909a_1   | conda-forge |
| libffi                 |   3.3             |     h58526e2_1   | conda-forge |
| libgcc                 |   7.2.0           |     h69d50b8_2   | conda-forge |
| libgcc-ng              |   9.3.0           |    h5dbcf3e_17   | conda-forge |
| libgfortran-ng         |   9.3.0           |    he4bcb1c_17   | conda-forge |
| libgfortran5           |   9.3.0           |    he4bcb1c_17   | conda-forge |
| libgomp                |   9.3.0           |    h5dbcf3e_17   | conda-forge |
| libnghttp2             |   1.41.0          |     h8cfc5f6_2   | conda-forge |
| libopenblas            |   0.3.12          | pthreads_h4812303_1   | conda-forge |
| libssh2                |   1.9.0           |     hab1572f_5   | conda-forge |
| libstdcxx-ng           |   9.3.0           |    h2ae2ef3_17   | conda-forge |
| ncurses                |   6.2             |     h58526e2_4   | conda-forge |
| openssl                |   1.1.1h          |     h516909a_0   | conda-forge |
| perl                   |   5.32.0          |     h36c2ea0_0   | conda-forge |
| pigz                   |   2.3.4           |     hed695b0_1   | conda-forge |
| pip                    |   20.3            |   pyhd8ed1ab_0   | conda-forge |
| python                 |   3.8.6           | hffdb5ce_1_cpython   | conda-forge |
| python_abi             |   3.8             |         1_cp38   | conda-forge |
| readline               |   8.0             |     he28a2e2_2   | conda-forge |
| samtools               |   1.7             |              1   | bioconda |
| setuptools             |   49.6.0          | py38h924ce5b_2   | conda-forge |
| sqlite                 |   3.34.0          |     h74cdb3f_0   | conda-forge |
| tbb                    |   2020.2          |     hc9558a2_0   | conda-forge |
| tk                     |   8.6.10          |     hed695b0_1   | conda-forge |
| wheel                  |   0.36.0          |   pyhd3deb0d_0   | conda-forge |
| xopen                  |   1.0.1           | py38h578d9bd_1   | conda-forge |
| xz                     |   5.2.5           |     h516909a_1   | conda-forge |
| zlib                   |   1.2.11          |  h516909a_1010   | conda-forge |




**Unfortunately, there's a bug in  the `Samtools` package of `conda` at the moment (i.e. `Samtools shared library libcrypto.so.1.0.0 not found`).**

**So we have to use the old version of `Samtools` in another `conda` enviroment (i.e., ` bifo-essential`), via**
```bash
conda activate bifo-essential
```
