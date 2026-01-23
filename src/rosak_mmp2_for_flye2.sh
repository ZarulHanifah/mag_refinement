#!/bin/bash

# Default values for parameters
PIDENT=95
COVERAGE=80
THREADS=1
END_REGION=2000
REFERENCE=""
READS=""
OUTPUT_PREFIX=""
TMPDIR=""

# Parse command-line arguments
while getopts "p:b:c:t:e:r:i:o:l:" opt; do
  case ${opt} in
    p) PIDENT=$OPTARG ;;
    b) TMPDIR=$OPTARG ;;
    c) COVERAGE=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    e) END_REGION=$OPTARG ;;
    r) REFERENCE=$OPTARG ;;
    i) READS=$OPTARG ;;
    o) OUTPUT_PREFIX=$OPTARG ;;
    l) LOGFILE=$OPTARG ;;
    *) echo "Usage: $0 [-p pident] [-c coverage] [-t threads] [-b tmpdir] [-e end_region] [-r reference.mmi] [-i reads.fastq.gz] [-o output_prefix] [-l logfile]"
       exit 1 ;;
  esac
done
shift $((OPTIND -1))

# Check if required parameters are provided
if [ -z "$REFERENCE" ] || [ ! -f "$REFERENCE" ]; then
    echo "Reference file $REFERENCE does not exist or is not specified." > "$LOGFILE"
    exit 1
fi

if [ -z "$READS" ] || [ ! -f "$READS" ]; then
    echo "Reads file $READS does not exist or is not specified." > "$LOGFILE"
    exit 1
fi

if [ -z "$TMPDIR" ]; then
    echo "TMPDIR $TMPDIR is not specified." > "$LOGFILE"
    exit 1
fi

if [ -z "$OUTPUT_PREFIX" ]; then
    echo "Output prefix is not specified." > "$LOGFILE"
    exit 1
fi

# Define output files
BAM_OUTPUT="${OUTPUT_PREFIX}_core.bam"
BAM_INDEX="${OUTPUT_PREFIX}_core.bam.bai"
IDS_OUTPUT="${OUTPUT_PREFIX}_core.ids"
FASTQ_OUTPUT="${OUTPUT_PREFIX}_core.fq"
END_BAM_OUTPUT="${OUTPUT_PREFIX}_ends.bam"
END_BAM_INDEX="${OUTPUT_PREFIX}_ends.bam.bai"
END_FASTQ_OUTPUT="${OUTPUT_PREFIX}_ends.fq"

# First alignment (Strict Filtering)
echo "##FIRST ALIGNMENT: High confidence alignments to contigs" > $LOGFILE
minimap2 -t "$THREADS" -ax map-ont "$REFERENCE" "$READS" 2>> $LOGFILE | \
samtools view -h -@ "$THREADS" | \
awk -v pident="$PIDENT" -v coverage="$COVERAGE" -v logfile="$LOGFILE" '
    BEGIN {OFS="\t"}
    /^@/ {print $0; next}  # Print headers
    {
        match($0, /NM:i:([0-9]+)/, arr);
        nm = arr[1];
        readID = $1 ;
        cigar = $6;
        read_len = length($10);

        match_len = 0;
        total_len = 0;

        while (match(cigar, /^[0-9]+[MIDNSHP=X]/)) {

            match(cigar, /([0-9]+)([MIDNSHP=X])/, matches);
            value = matches[1] + 0;
            op = matches[2];

            if (op == "M") {
                match_len += value;
            }
            if (op != "H" && op != "S") {
                total_len += value;
            }

           cigar = substr(cigar, length(matches[0]) + 1);
        }

        if (match_len == 0 || read_len <= 0 || !nm || total_len == 0) {
            # print "Discarded (invalid alignment):", readID >> logfile;
            next;
        }

        # Ensure match_len does not exceed read_len
        match_len = (match_len > read_len) ? read_len : match_len;

        #percent_identity = (match_len - nm) / match_len * 100;
        #alignment_coverage = match_len / read_len * 100;

        percent_identity = (match_len - nm) / total_len * 100;
        alignment_coverage = total_len / read_len * 100;

        if (percent_identity < pident || alignment_coverage < coverage) {
            # print "Discarded:" readID, "readLen=" read_len, "cigar=" $6, "pident=" percent_identity, "coverage=" alignment_coverage >> logfile;
            next;
        }
        else {
            print "Pass:", readID, "readLen=" read_len, "cigar=" $6, "pident=" percent_identity, "coverage=" alignment_coverage >> logfile;
            print $0;
            next;
        }
    }' | \
samtools sort -@ "$THREADS" -T $TMPDIR 2>> $LOGFILE | \
samtools view | \ 
cut -f1 | sort -u > $IDS_OUTPUT

seqtk subseq $READS $IDS_OUTPUT > $FASTQ_OUTPUT

#samtools index -@ "$THREADS" "$BAM_OUTPUT" "$BAM_INDEX"
#samtools bam2fq -@ "$THREADS" "$BAM_OUTPUT" > "$FASTQ_OUTPUT" 2>> $LOGFILE


# Second alignment (Lenient Filtering for Contig Ends)
#echo ; echo "##SECOND ALIGNMENT: Align toward contig ends" >> $LOGFILE
#minimap2 -t "$THREADS" -ax map-ont "$REFERENCE" "$READS" 2>> $LOGFILE | \
#samtools view -h -@ "$THREADS" | \
#awk -v end_region="$END_REGION" -v logfile="$LOGFILE" -v pident="$PIDENT" '
#    BEGIN {
#        OFS = "\t";
#    }
#
#    # Store contig lengths from the SAM header (@SQ lines)
#    /^@SQ/ {
#        split($2, sn_field, ":");  # Extract SN (contig name)
#        split($3, ln_field, ":");  # Extract LN (contig length)
#        contig_lengths[sn_field[2]] = ln_field[2];  # Store in associative array
#        print $0;
#        next;
#    }
#
#    # Print header lines without modification
#    /^@/ {
#        print $0;
#        next;
#    }
#
#    # Process alignment records
#    {
#        contig_name = $3;               # Contig name (reference sequence)
#        start_pos = $4;                 # Read alignment start position
#        read_length = length($10);      # Length of the query sequence
#
#        # Get the contig length from the stored header values
#        contig_length = contig_lengths[contig_name];
#
#        # Extract NM (edit distance) field from optional fields
#        nm = -1;
#        for (i = 12; i <= NF; i++) {
#            if ($i ~ /^NM:i:/) {
#                split($i, nm_field, ":");
#                nm = nm_field[3];  # NM value (edit distance)
#                break;
#            }
#        }
#
#        # Calculate percent identity
#        if (nm >= 0) {
#            matches = read_length - nm;  # Matches = alignment length - mismatches
#            percent_identity = (matches / read_length) * 100;
#        } else {
#            percent_identity = 0;  # Default to 0 if NM field is missing
#        }
#
#
#
#
#        # If contig length is found, apply the filtering condition
#        if (contig_length && 
#            (read_length >= 2000) && 
#            !(start_pos <= 2000 && (start_pos + read_length) >= contig_length - 2000) &&
#            (start_pos >= 2000 || (start_pos + read_length <= contig_length - 2000)) &&
#            (percent_identity >= pident)) {
#            print $0;
#        }
#    }' | \
#samtools sort -@ "$THREADS" -T $TMPDIR | \
#samtools view -b -@ "$THREADS" -o "$END_BAM_OUTPUT"
#
#samtools index -@ "$THREADS" "$END_BAM_OUTPUT" "$END_BAM_INDEX"
#samtools bam2fq -@ "$THREADS" "$END_BAM_OUTPUT" > "$END_FASTQ_OUTPUT" 2>> $LOGFILE
#
#echo "End-region BAM saved to $END_BAM_OUTPUT" >> $LOGFILE
#echo "BAM index saved to $END_BAM_INDEX" >> $LOGFILE
#echo "End-region FASTQ saved to $END_FASTQ_OUTPUT" >> $LOGFILE

