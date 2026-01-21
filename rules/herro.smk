rule minimap2_index:
    input: config["dereplicated_genome_path"]
    output: temp(os.path.join(temp_path, "minimap2_index/{mag}/{mag}.mmi"))
    log: os.path.join(results_path, "log/minimap2_index/{mag}.log")
    threads: 8
    shell:
        """ 
        module load minimap2/2.28
        minimap2 -t {threads} -x map-ont -d {output} {input}
        """ 

def get_osample_input_directory(wildcards):
    return samples_dorado[wildcards.osample]

checkpoint chunk_osample:
    output: os.path.join(results_path, "minimap2/chunking/{osample}.map")
    threads: 1
    params:
        input_folder = get_osample_input_directory,
    shell:
        """ 
        python src/chunking_fastq.py -i {params.input_folder} | sed "s/^chunk/{wildcards.osample}_chunk/" > {output}
        """ 

def output_fastq_chunk_osample(wildcards):
    return expand(checkpoints.chunk_osample.get(**wildcards).output, osample=samples_list)

rule get_fastq_chunk:
    input:
        chunkmap = output_fastq_chunk_osample
    output:
        temp(os.path.join(temp_path, "minimap2/fastq_depot/{osample}/{fastq}.fastq.gz"))
    shell:
        """ 
        cat $(grep -w {wildcards.fastq} {input.chunkmap} | cut -f2 | tr "," " ") > {output}
        """ 

# Convert minimap2_mini from a checkpoint to a rule
rule minimap2_mini:
    input:
        reads = rules.get_fastq_chunk.output,
        idx = rules.minimap2_index.output
    output:
        bam = temp(os.path.join(temp_path, "minimap2_mini/{mag}/{osample}/{fastq}.bam")),
        bai = temp(os.path.join(temp_path, "minimap2_mini/{mag}/{osample}/{fastq}.bam.bai")),
        fq = temp(os.path.join(temp_path, "minimap2_mini/{mag}/{osample}/{fastq}.fq"))
    log: os.path.join(results_path, "log/minimap2_mini/{mag}/{osample}/{fastq}.log")
    threads: 2
    params:
        pident = config["pident"],
        coverage = config["coverage"]
    shell:
        """ 
        module load minimap2/2.28
        module load samtools
    
        PREFIX=$(echo {output.fq} | sed "s/\.fq//")

        tmpdir=$(dirname {output.fq})"_"{wildcards.fastq}"/"
        mkdir -p $tmpdir

        ./src/minimap2_stringent.sh -t {threads} -p {params.pident} -c {params.coverage} -r {input.idx} -i {input.reads} -o $PREFIX -b $tmpdir -l {log}

        rm -rf $tmpdir
        """ 

rule minimap2_fq:
    input:
        reads = rules.get_fastq_chunk.output,
        idx = rules.minimap2_index.output
    output:
        fq = os.path.join(temp_path, "minimap2_fq/{mag}/{osample}/{fastq}.fq"),
        ids = temp(os.path.join(temp_path, "minimap2_fq/{mag}/{osample}/{fastq}.ids"))
    log: os.path.join(results_path, "log/minimap2_fq/{mag}/{osample}/{fastq}.log")
    conda: "seqkit_"
    threads: 2
    params:
        pident = 95,
        coverage = 50,
        end_len = 2000
    shell:
        """ 
        module load minimap2/2.28
        module load samtools
        module load seqtk
    
        PREFIX=$(echo {output.fq} | sed "s/\.fq//")

        tmpdir=$(dirname {output.fq})"_"{wildcards.fastq}"/"
        mkdir -p $tmpdir

         # -e {params.end_len} \

        ./src/mmp2_for_flye.sh -p {params.pident} \
         -c {params.coverage} \
         -t {threads} \
         -r {input.idx} \
         -i {input.reads} \
         -o $PREFIX \
         -b $tmpdir -l {log}
        
        # cat $PREFIX"_core.fq" $PREFIX"_ends.fq" | seqkit rmdup  > {output.fq} 
        # cat $PREFIX"_core.fq" | seqkit rmdup  > {output.fq} 

        rm -rf $tmpdir
        """ 

# rule herro:
#     input:
#         fq = rules.concat_fastq_aligned_to_mag.output,
#         herro_model = config["herro_model"]
#     output:
#         preprocess = os.path.join(results_path, "herro/preprocess/{mag}.fastq.gz"),
#         ids = os.path.join(results_path, "herro/preprocess/{mag}.ids"),
#         batch_alignment = directory(os.path.join(results_path, "herro/batch_alignment/{mag}")),
#         fasta= os.path.join(results_path, "herro/fasta/{mag}.fasta"),
#     conda: "seqkit_"
#     threads: 8
#     params:
#         gpu = 0,
#         batch_size = 32
#     shell:
#         """
#         preprocess_prefix=$(echo {output.preprocess} | sed "s/\.fastq.gz//")
#         module load herro
#         module load minimap2/2.24
#         
#         rm -rf $(dirname {output.preprocess})
#         src/preprocess.sh {input.fq} $preprocess_prefix {threads} 1
# 
#         cat {output.preprocess} | seqkit seq -n -i > {output.ids}
#         src/create_batched_alignments.sh {output.preprocess} {output.ids} {threads} {output.batch_alignment}
# 
#         herro inference --read-alns {output.batch_alignment} -t {threads} \
#          -d {params.gpu} -m {input.herro_model} -b {params.batch_size} \
#          {output.preprocess} {output.fasta}
#         """
