#!/usr/bin/env python

import argparse
import sys


def parse_args():
    p = argparse.ArgumentParser(description="Convert GFA to FASTA (extract S lines)")
    p.add_argument("--gfa", help="Input GFA file")
    p.add_argument("-o", "--output", default="-", help="Output FASTA file (default: stdout)")
    return p.parse_args()

def main():
    args = parse_args()

    out = sys.stdout if args.output == "-" else open(args.output, "w")

    with open(args.gfa) as f:
        for line in f:
            if not line.startswith("S\t"):
                continue
            fields = line.strip().split("\t")
            name = fields[1]
            seq = fields[2]
            if seq == "*":
                continue
        out.write(f">{name}\n")
        out.write(seq + "\n")

    if out is not sys.stdout:
        out.close()

if __name__ == "__main__":
    main()
