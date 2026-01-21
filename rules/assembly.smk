def get_samples_for_mag(wildcards):
    mag_name = wildcards.mag
    try:
        mag = sesh.get_mag(mag_name)
        abund_dict = mag.get_depth_df().sum().to_dict()
        present_samples = [
            sample_name
            for sample_name, abundance in abund_dict.items()
            if abundance > 0.5
        ]
        return present_samples
    except FileNotFoundError as e:
        raise e(f"Warning: MAG file for {mag_name} not found. Skipping.")

def get_fastq_chunks(osample):
    """ 
    Fetches chunk identifiers for fastq files from the chunk_osample checkpoint output.
    """ 
    chunk_map_file = str(checkpoints.chunk_osample.get(osample=osample).output)
    df = pd.read_csv(chunk_map_file, sep="\t", header=None, index_col=0)
    return df.index.tolist()

def get_all_fq(wildcards):
    """ 
    Collect all fastq outputs from the minimap2_mini rule for a given mag.
    """ 
    fastqs = []

    relevant_samples = get_samples_for_mag(wildcards)

    for osample in relevant_samples:
        fastq_chunks = get_fastq_chunks(osample)
        for fastq in fastq_chunks:
            # Dynamically retrieve outputs from the rule
            fq_output = os.path.join(temp_path, "minimap2_fq/{mag}/{osample}/{fastq}.fq").format(mag=wildcards.mag, osample=osample, fastq=fastq)
            fastqs.append(fq_output)
    return fastqs

# temp(directory(os.path.join(temp_path, "minimap2_fq_collected/{mag}")))

rule collect_minimap2_fq:
    input:
        lambda wildcards: get_all_fq(wildcards)
    output:
        temp(os.path.join(temp_path, "minimap2_fq_collected/{mag}/all.fq"))
    shell:
        """
        # mkdir -p {output}
        # cat {input} > {output}/all.fq

        cat {input} > {output}
        """

rule flye_fq:
    input:
        fq = lambda wc: os.path.join(temp_path, f"minimap2_fq_collected/{wc.mag}/all.fq")
    output:
        assem = os.path.join(results_path, "flye_fq/{mag}/assembly.fasta")
    log:
        os.path.join(results_path, "log/flye_fq/{mag}.log")
    conda: "flye_"
    threads: 8
    shell:
        """
        outdir=$(dirname {output.assem})

        flye --nano-hq {input.fq} --out-dir $outdir --threads {threads} --meta 2> {log}
        """

rule hifiasm_fq:
    input:
        fq = lambda wc: os.path.join(temp_path, f"minimap2_fq_collected/{wc.mag}/all.fq")
    output:
        assem = os.path.join(results_path, "hifiasm_fq/{mag}/{mag}.p_ctg.gfa")
    log:
        os.path.join(results_path, "log/hifiasm_fq/{mag}.log")
    conda: "hifiasm_"
    threads: 4
    shell:
        """
        prefix=$(echo {output.assem}| sed "s/\.p_ctg\.gfa//")

        hifiasm --ont -o $prefix -t {threads} --primary --n-hap 1  {input.fq} 2> {log}
        """

rule myloasm_fq:
    input:
        fq = lambda wc: os.path.join(temp_path, f"minimap2_fq_collected/{wc.mag}/all.fq")
    output:
        assem = os.path.join(results_path, "myloasm_fq/{mag}/assembly_primary.fa")
    conda: "myloasm_"
    log:
        os.path.join(results_path, "log/myloasm_fq/{mag}.log")
    threads: 8
    shell:
        """
        outdir=$(dirname {output.assem})
        rm -rf $outdir

        myloasm {input.fq} --output-dir $outdir --threads {threads} --clean-dir &> {log}
        """

rule wtdbg2_fq:
    input:
        fq = lambda wc: os.path.join(temp_path, f"minimap2_fq_collected/{wc.mag}/all.fq")
    output:
        assem = os.path.join(results_path, "wtdbg2_fq/{mag}/{mag}.ctg.gfa")
    log:
        os.path.join(results_path, "log/wtdbg2_fq/{mag}.log")
    threads: 4
    shell:
        """
        module load wtdbg2
        prefix=$(echo {output.assem}| sed "s/\.ctg\.gfa//")
        
        wtdbg2 -t {threads} -i {input.fq} -fo $prefix 2> {log}
        wtpoa-cns -t {threads} -i $prefix".ctg.lay.gz" -fo $prefix".raw.fa" 2> {log}
        ./src/wtdbg-dot2gfa.pl $prefix".ctg.dot.gz" {output} 2> {log}
        """
