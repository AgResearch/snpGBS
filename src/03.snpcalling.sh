#!/bin/bash
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
