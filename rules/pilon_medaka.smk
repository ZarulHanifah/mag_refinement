rule minimap2_index_pilon_medaka:
    input: rules.medaka.output.assem
    output: temp(os.path.join(temp_path, "pilon/medaka/minimap2_index/{mag}/{mag}.mmi"))
    log: os.path.join(results_path, "log/pilon/medaka/minimap2_index/{mag}.log")
    threads: 8
    shell:
        """ 
        module load minimap2/2.28
        minimap2 -t {threads} -x map-ont -d {output} {input}
        """ 

rule align_short_read_per_sample_pilon_medaka:
    input:
        idx = rules.minimap2_index_pilon_medaka.output,
        fwd = lambda wildcards: config['extra_short_reads'][wildcards.sr_sample]['fwd'],
        rvs = lambda wildcards: config['extra_short_reads'][wildcards.sr_sample]['rvs'],
    output:
        bam = temp(os.path.join(temp_path, 'pilon/medaka/sr_alignment/{mag}/{sr_sample}.bam')),
        bai = temp(os.path.join(temp_path, 'pilon/medaka/sr_alignment/{mag}/{sr_sample}.bam.bai')),
        flagstat = os.path.join(temp_path, 'pilon/medaka/sr_alignment/{mag}/{sr_sample}.flagstat.txt'),
        fq1 = temp(os.path.join(temp_path, 'pilon/medaka/sr_alignment/{mag}/{sr_sample}_R1.fastq.gz')),
        fq2 = temp(os.path.join(temp_path, 'pilon/medaka/sr_alignment/{mag}/{sr_sample}_R2.fastq.gz')),
        fqu = temp(os.path.join(temp_path, 'pilon/medaka/sr_alignment/{mag}/{sr_sample}_unpaired.fastq.gz')),
    log:
        os.path.join(results_path, 'log/pilon/medaka/align_short_read_per_sample_pilon_medaka/{mag}__{sr_sample}.log')
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

rule merge_and_split_bams_for_pilon_medaka:
    input:
        bams = lambda wildcards: expand(rules.align_short_read_per_sample_pilon_medaka.output.bam, mag=wildcards.mag, sr_sample=sr_samples_list),
        bais = lambda wildcards: expand(rules.align_short_read_per_sample_pilon_medaka.output.bai, mag=wildcards.mag, sr_sample=sr_samples_list)
    output:
        p_bam = temp(os.path.join(temp_path, 'pilon/medaka/paired/{mag}.bam')),
        p_bai = temp(os.path.join(temp_path, 'pilon/medaka/paired/{mag}.bam.bai')),
        up_bam = temp(os.path.join(temp_path, 'pilon/medaka/unpaired/{mag}.bam')),
        up_bai = temp(os.path.join(temp_path, 'pilon/medaka/unpaired/{mag}.bam.bai')),
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

rule pilon_medaka:
    input:
        mag = rules.medaka.output.assem,
        p_bam = rules.merge_and_split_bams_for_pilon_medaka.output.p_bam,
        p_bai = rules.merge_and_split_bams_for_pilon_medaka.output.p_bai,
        up_bam = rules.merge_and_split_bams_for_pilon_medaka.output.up_bam,
        up_bai = rules.merge_and_split_bams_for_pilon_medaka.output.up_bai,
    output:
        assem = os.path.join(results_path, 'pilon/medaka/{mag}/{mag}.fasta'),
    log: os.path.join(results_path, 'log/pilon/medaka/{mag}.log')
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


