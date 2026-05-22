rule align_short_read_per_sample:
    input:
        idx = rules.minimap2_index.output,
        fwd = lambda wildcards: config['extra_short_reads'][wildcards.sr_sample]['fwd'],
        rvs = lambda wildcards: config['extra_short_reads'][wildcards.sr_sample]['rvs'],
    output:
        bam = temp(os.path.join(temp_path, 'racon/sr_alignment/{mag}/{sr_sample}.bam')),
        bai = temp(os.path.join(temp_path, 'racon/sr_alignment/{mag}/{sr_sample}.bam.bai')),
        fq1 = temp(os.path.join(temp_path, 'racon/sr_alignment/{mag}/{sr_sample}_R1.fastq.gz')),
        fq2 = temp(os.path.join(temp_path, 'racon/sr_alignment/{mag}/{sr_sample}_R2.fastq.gz')),
    log:
        os.path.join(results_path, 'log/racon/align_short_read_per_sample/{mag}__{sr_sample}.log')
    threads: 4
    shell:
        '''
        module load minimap2/2.28
        module load samtools/1.19.3

        tmpdir=$(dirname {output.bam})_tmp
        outdir=$(dirname {output.bam})
        rm -rf $tmpdir ; mkdir -p $tmpdir

        minimap2 -ax sr -t {threads} {input.idx} {input.fwd} {input.rvs} 2> {log} | \
         samtools view -@ {threads} -h -F 4 - | \
         samtools sort -@ {threads} -T $tmpdir -o {output.bam} - 2>> {log}
        samtools index {output.bam}

        samtools fastq -@ {threads} -1 {output.fq1} -2 {output.fq2} -0 /dev/null -s /dev/null {output.bam} 2>> {log}

        rm -rf $tmpdir
        '''

rule merge_sr_alignment:
    input:
        bams = lambda wildcards: expand(rules.align_short_read_per_sample.output.bam, mag=wildcards.mag, sr_sample=sr_samples_list),
        bais = lambda wildcards: expand(rules.align_short_read_per_sample.output.bai, mag=wildcards.mag, sr_sample=sr_samples_list)
    output:
        temp(os.path.join(temp_path, 'racon/merged_alignment/{mag}.sam'))
    threads: 4
    shell:
        '''
        module load minimap2/2.28
        module load samtools/1.19.3
        samtools merge -@ {threads} -O SAM {output} {input.bams}
        '''

rule racon:
    input:
        mag = ancient(config["dereplicated_genome_path"]),
        sam = rules.merge_sr_alignment.output,
        fq1 = lambda wildcards: expand(rules.align_short_read_per_sample.output.fq1, mag=wildcards.mag, sr_sample=sr_samples_list),
        fq2 = lambda wildcards: expand(rules.align_short_read_per_sample.output.fq2, mag=wildcards.mag, sr_sample=sr_samples_list)
    output:
        fq = temp(os.path.join(temp_path, 'racon/merged_alignment/{mag}.fastq.gz')),
        assem = os.path.join(results_path, 'racon/{mag}.fasta')
    threads: 8
    log:
        os.path.join(results_path, 'log/racon/racon/{mag}.log')
    shell:
        '''
        module load racon/1.5.0

        zcat {input.fq1} {input.fq2} | gzip > {output.fq}

        racon -t {threads} -m 8 -x -6 -g -8 -w 500 \
         {output.fq} \
         {input.sam} \
         {input.mag} > {output.assem} 2> {log}
        '''
