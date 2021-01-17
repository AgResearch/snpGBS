#!/bin/bash
filename='barcodes.txt'
n=0
while read line; do
 echo ">barcode$n"
 echo "^$line"
 n=$((n+1))
done < $filename
