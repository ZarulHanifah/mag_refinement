#!/usr/bin/bash
#SBATCH --job-name=ZAR_debug         # Job name
#SBATCH --cpus-per-task=4            # Number of CPU cores per task
#SBATCH --output=my_job.out          # Standard output file
#SBATCH --error=my_job.err           # Standard error file
#SBATCH --time=01:00:00              # Time limit (hh:mm:ss)
#SBATCH --mem=8G                     # Memory per node
#SBATCH --partition=highprio

module load minimap2
module load samtools
module load seqtk

PREFIX=$(echo results/minimap2_fq/C3D5_F_metabat.1225/C3D5/C3D5_chunk10.fq | sed "s/\.fq//")

outdir="results/minimap2_fq/C3D5_F_metabat.1225/C3D5/"

tmpdir=$(dirname results/minimap2_fq/C3D5_F_metabat.1225/C3D5/C3D5_chunk10_core.bam)"_"tmp"/"
rm -rf $outdir ; mkdir -p $outdir
mkdir -p $tmpdir

#./src/minimap2_stringent.sh -p 95 -c 50 -t 4 -r input_folder/dereplicated_genomes/C3D5_F_metabat.1225.mmi -i results/minimap2/fastq_depot/C3D5/C3D5_chunk10.fastq.gz -o $PREFIX -b $tmpdir -l results/log/minimap2_fq/C3D5_F_metabat.1225/C3D5/C3D5_chunk10.log

./src/mmp2_for_flye.sh -p 95 -c 50 -t 4 -r input_folder/dereplicated_genomes/C3D5_F_metabat.1225.mmi -i results/minimap2/fastq_depot/C3D5/C3D5_chunk10.fastq.gz -o $PREFIX -b $tmpdir -l results/log/minimap2_fq/C3D5_F_metabat.1225/C3D5/C3D5_chunk10.log

#./src/mmp2_for_flye2.sh -p 95 -c 50 -t 4 -e 2000 -r input_folder/dereplicated_genomes/C3D5_F_metabat.1225.mmi -i results/minimap2/fastq_depot/C3D5/C3D5_chunk10.fastq.gz -o $PREFIX -b $tmpdir -l results/log/minimap2_fq/C3D5_F_metabat.1225/C3D5/C3D5_chunk10.log

# cat $PREFIX"_core.fq" $PREFIX"_ends.fq" > results/minimap2_fq/C3D5_F_metabat.1225/C3D5/C3D5_chunk10.fq

rm -rf $tmpdir

