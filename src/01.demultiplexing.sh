#!/bin/bash
cutadapt -j 8 -e 0.15 --no-indels -g file:barcodes.fasta -o "demultiplexed-{name}.fastq.gz" example.fastq >01.demultiplexed.stdout 2>01.demultiplexed.stderr

## -j: number of CPU cores to be used
## -e: maximum error rate to 15%
## --no-indels: insertions or deletions are not allowed
## -g:  barcode adapters are attcaed to the 5’ end
## file:barcodes.fasta: file contains a list of barcodes
## -o: output file name pattern
## input file at the end (e.g., example.fastq)
