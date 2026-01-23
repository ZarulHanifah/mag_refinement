#!/bin/bash
# Convert a multi-line FASTA file into a single line FASTA file
# These are easier/faster to process using native UNIX tools like paste - - < in.fasta
 awk 'BEGIN { RS = "\n>"; FS = "\n"; OFS = "" };
{
  if (NR == 1) {
    print $1
  } else 
  if (NR > 1) {
    print ">"$1
  }
  $1=""
  print
}' 
