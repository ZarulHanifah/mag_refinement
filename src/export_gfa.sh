#!/usr/bin/env bash

# Parse command-line arguments
while getopts "i:" opt; do
  case ${opt} in
    i) input_folder=$OPTARG ;;
    *) echo "Usage: $0 [-i input_folder]"
       exit 1 ;;
  esac
done

# Ensure the input_folder is specified
if [ -z "$input_folder" ]; then
    echo "Input folder is not specified. Use -i to provide it."
    exit 1
fi

# Loop through files and generate rsync commands
for i in $(find "$input_folder" -type f | grep "assembly_graph.gfa"); do
    mag=$(basename $(dirname $i))
    path=$(readlink -f $i)
    echo "rsync -r --progress mzar0002@montage:$path $mag.gfa"
done

