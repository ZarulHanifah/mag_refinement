#!/bin/bash

# Usage: ./mmp2_for_flye.sh -p 95 -c 50 -t 4 -r reference.mmi -i reads.fastq.gz -o output_prefix -b tmpdir -l logfile

# Default values for parameters
PIDENT=95
COVERAGE=80
THREADS=1
END_REGION=2000
REFERENCE=""
READS=""
OUTPUT_PREFIX=""
TMPDIR=""
LOGFILE=""

# Parse command-line arguments
while getopts "p:b:c:e:t:r:i:o:l:" opt; do
  case ${opt} in
    p) PIDENT=$OPTARG ;;
    b) TMPDIR=$OPTARG ;;
    c) COVERAGE=$OPTARG ;;
    e) END_REGION=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    r) REFERENCE=$OPTARG ;;
    i) READS=$OPTARG ;;
    o) OUTPUT_PREFIX=$OPTARG ;;
    l) LOGFILE=$OPTARG ;;
    *) echo "Usage: $0 [-p pident] [-b tmpdir] [-c coverage] [-e end_region]  [-t threads] [-r reference.mmi] [-i reads.fastq.gz] [-o output_prefix] [-l logfile]"
       exit 1 ;;
  esac
done
shift $((OPTIND -1))

# Check if reference, reads, output prefix, and logfile are provided
if [ -z "$REFERENCE" ] || [ ! -f "$REFERENCE" ]; then
    echo "Reference file $REFERENCE does not exist or is not specified." >> "$LOGFILE"
    exit 1
fi

if [ -z "$READS" ] || [ ! -f "$READS" ]; then
    echo "Reads file $READS does not exist or is not specified." >> "$LOGFILE"
    exit 1
fi

if [ -z "$TMPDIR" ]; then
    echo "TMPDIR $TMPDIR is not specified." >> "$LOGFILE"
    exit 1
fi

if [ -z "$OUTPUT_PREFIX" ]; then
    echo "Output prefix is not specified." >> "$LOGFILE"
    exit 1
fi

if [ -z "$LOGFILE" ]; then
    echo "Log file is not specified." >> "$LOGFILE"
    exit 1
fi

# Define output files
BAM_OUTPUT="${OUTPUT_PREFIX}.bam"
BAM_INDEX="${OUTPUT_PREFIX}.bam.bai"
IDS_OUTPUT="${OUTPUT_PREFIX}.ids"
FASTQ_OUTPUT="${OUTPUT_PREFIX}.fq"

# Run minimap2 using specified number of threads and pipe through samtools
echo "##FIRST ALIGNMENT: High confidence alignments to contigs" > $LOGFILE
minimap2 -t "$THREADS" -ax map-ont "$REFERENCE" "$READS" | \
samtools view -h -@ "$THREADS" | \
awk -v pident="$PIDENT" -v coverage="$COVERAGE" -v logfile="$LOGFILE" -v end_region="$END_REGION" '
    BEGIN {OFS="\t"}
    /^@/ {print $0; next}  # Print header lines
    {
        nm = 0;
        match($0, /NM:i:([0-9]+)/, arr);
        if (arr[1] != "") nm = arr[1] + 0;
        cigar = $6;

        # Calculate match length and total alignment length
        match_len = 0;
        total_len = 0;  # Total alignment length

        while (match(cigar, /^[0-9]+[MIDNSHP=X]/)) {
            match(cigar, /([0-9]+)([MIDNSHP=X])/, matches);
            value = matches[1] + 0;  # Convert to integer
            op = matches[2];         # Operation type

            if (op == "M") {
                match_len += value;
            }
            # if (op != "H" && op != "S") {
            #     total_len += value;
            # }
            total_len += value;

            # cigar = substr(cigar, length(matches[0]) + 1);
            cigar = substr(cigar, RLENGTH + 1);
        }

        read_len = length($10);  # Sequence length

        # if (match_len == 0 || read_len <= 0 || !nm || total_len == 0) {
	
	# Removed !nm to accept perfect alignments
        if (match_len == 0 || read_len <= 0 || total_len == 0) {
            next;
        }

        # Ensure match_len does not exceed read_len
        match_len = (match_len > read_len) ? read_len : match_len;

        percent_identity = (match_len - nm) / match_len * 100;
        alignment_coverage = match_len / read_len * 100;

        if (read_len < end_region || percent_identity < pident || alignment_coverage < coverage) {
            # print "Discarded:", $1, "pident=" percent_identity, "coverage=" alignment_coverage >> logfile;
            next;
        } else {
            print "Pass:", $1, "readLen=" read_len, "pident=" percent_identity, "coverage=" alignment_coverage, "cigar=" $6 >> logfile;
            print $0;
        }
    }'  |
samtools sort -@ $THREADS -T $TMPDIR 2>> $LOGFILE | \
tee $BAM_OUTPUT | \
samtools view -@ $THREADS | cut -f1 | sort -u > $IDS_OUTPUT

samtools index -@ $THREADS $BAM_OUTPUT $BAM_INDEX

seqtk subseq $READS $IDS_OUTPUT > $FASTQ_OUTPUT

echo "Filtered BAM saved to $BAM_OUTPUT" >> $LOGFILE
echo "BAM index saved to $BAM_INDEX" >> $LOGFILE
echo "Filtered FASTQ saved to $FASTQ_OUTPUT" >> $LOGFILE
