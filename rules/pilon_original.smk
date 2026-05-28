rule minimap2_index_pilon_original:
    input: ancient(config["dereplicated_genome_path"])
    output: temp(os.path.join(temp_path, "pilon/original/minimap2_index/{mag}/{mag}.mmi"))
    log: os.path.join(results_path, "log/pilon/original/minimap2_index/{mag}.log")
    threads: 8
    shell:
        """ 
        module load minimap2/2.28
        minimap2 -t {threads} -x map-ont -d {output} {input}
        """ 

rule align_short_read_per_sample_pilon_original:
    input:
        idx = rules.minimap2_index_pilon_original.output,
        fwd = lambda wildcards: config['extra_short_reads'][wildcards.sr_sample]['fwd'],
        rvs = lambda wildcards: config['extra_short_reads'][wildcards.sr_sample]['rvs'],
    output:
        bam = temp(os.path.join(temp_path, 'pilon/original/sr_alignment/{mag}/{sr_sample}.bam')),
        bai = temp(os.path.join(temp_path, 'pilon/original/sr_alignment/{mag}/{sr_sample}.bam.bai')),
        flagstat = os.path.join(temp_path, 'pilon/original/sr_alignment/{mag}/{sr_sample}.flagstat.txt'),
        fq1 = temp(os.path.join(temp_path, 'pilon/original/sr_alignment/{mag}/{sr_sample}_R1.fastq.gz')),
        fq2 = temp(os.path.join(temp_path, 'pilon/original/sr_alignment/{mag}/{sr_sample}_R2.fastq.gz')),
        fqu = temp(os.path.join(temp_path, 'pilon/original/sr_alignment/{mag}/{sr_sample}_unpaired.fastq.gz')),
    log:
        os.path.join(results_path, 'log/pilon/original/align_short_read_per_sample_pilon_original/{mag}__{sr_sample}.log')
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
        samtools flagstat {output.bam} > {output.flagstat}

        samtools fastq -@ {threads} -1 {output.fq1} -2 {output.fq2} -0 {output.fqu} -s /dev/null {output.bam} 2>> {log}

        rm -rf $tmpdir
        '''

rule merge_and_split_bams_for_pilon_original:
    input:
        bams = lambda wildcards: expand(rules.align_short_read_per_sample_pilon_original.output.bam, mag=wildcards.mag, sr_sample=sr_samples_list),
        bais = lambda wildcards: expand(rules.align_short_read_per_sample_pilon_original.output.bai, mag=wildcards.mag, sr_sample=sr_samples_list)
    output:
        p_bam = temp(os.path.join(temp_path, 'pilon/original/paired/{mag}.bam')),
        p_bai = temp(os.path.join(temp_path, 'pilon/original/paired/{mag}.bam.bai')),
        up_bam = temp(os.path.join(temp_path, 'pilon/original/unpaired/{mag}.bam')),
        up_bai = temp(os.path.join(temp_path, 'pilon/original/unpaired/{mag}.bam.bai')),
    threads: 4
    shell:
        '''
        module load minimap2/2.28
        module load samtools/1.19.3

        samtools merge -@ {threads} -o - {input.bams} | \
         samtools view -@ {threads} -b -f 1 -o {output.p_bam} -
        samtools index {output.p_bam}

        samtools merge -@ {threads} -o - {input.bams} | \
         samtools view -@ {threads} -b -G 1 -o {output.up_bam} -
        samtools index {output.up_bam}
        '''

rule pilon_original:
    input:
        mag = ancient(config["dereplicated_genome_path"]),
        p_bam = rules.merge_and_split_bams_for_pilon_original.output.p_bam,
        p_bai = rules.merge_and_split_bams_for_pilon_original.output.p_bai,
        up_bam = rules.merge_and_split_bams_for_pilon_original.output.up_bam,
        up_bai = rules.merge_and_split_bams_for_pilon_original.output.up_bai,
    output:
        assem = os.path.join(results_path, 'pilon/original/{mag}/{mag}.fasta'),
    log: os.path.join(results_path, 'log/pilon/original/{mag}.log')
    threads: 4
    shell:
        '''
        module load pilon/1.22

        outdir=$(dirname {output.assem})

        java -Xmx60g -jar /usr/local/pilon/1.22/pilon-1.22.jar \
         --threads {threads} \
         --genome {input.mag} \
         --frags {input.p_bam} \
         --unpaired {input.up_bam} \
         --outdir $outdir \
         --output {wildcards.mag} \
         --verbose --changes --tracks &> {log}
        '''

