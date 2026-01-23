#!/bin/bash

# Loop over all mag folders
for mag_folder in results/magpurify2_coverage/*/; do
    # Check if the folder contains 'coverage_scores.tsv'
    if [[ ! -f "$mag_folder/scores/coverage_scores.tsv" ]]; then
        # If the file does not exist, remove the mag folder
        echo "Removing $mag_folder"
        rm -r "$mag_folder"
    fi
done
