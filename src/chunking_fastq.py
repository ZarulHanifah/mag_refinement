import os
import subprocess
import argparse

# Function to get the size of a file
def get_file_size(file_path):
    return os.path.getsize(file_path)

def main(input_folder, chunk_size_gb, output_file=None):
    size_threshold = chunk_size_gb * 1024**3  # Convert chunk size from GB to bytes

    # Get the list of files
    files = subprocess.check_output(f'find $(readlink -f {input_folder}) -type f | grep "fastq.gz"', shell=True).decode().splitlines()

    # Calculate the size of each file
    file_sizes = [(file, get_file_size(file)) for file in files]

    # Sort files by size (optional, but can help optimize the concatenation)
    file_sizes.sort(key=lambda x: x[1], reverse=True)

    chunks = []
    current_chunk = []
    current_size = 0

    for file, size in file_sizes:
        if current_size + size > size_threshold and current_chunk:
            chunks.append(current_chunk)
            current_chunk = []
            current_size = 0

        current_chunk.append(file)
        current_size += size

    if current_chunk:
        chunks.append(current_chunk)

    # Output the result
    output = ""
    for i, chunk in enumerate(chunks):
        output += f"chunk{i+1}\t{','.join(chunk)}\n"

    if output_file:
        with open(output_file, "w") as f:
            f.write(output)
        print(f"Output table generated as '{output_file}'.")
    else:
        print(output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate files into chunks of at least the specified size.")
    parser.add_argument("-i", "--input_folder", required=True, help="Input folder containing the files.")
    parser.add_argument("-s", "--size", type=int, default=5, help="Desired chunk size in GB. Default is 5 GB.")
    parser.add_argument("-o", "--output_file", help="Output file to save the result. If not provided, output will be printed to stdout.")

    args = parser.parse_args()
    main(args.input_folder, args.size, args.output_file)

