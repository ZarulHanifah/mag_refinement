def get_bam_osample(wildcards):
    fastqs = get_fastq_chunks(wildcards.osample)
    return expand(rules.minimap2_fq.output.bam, mag=wildcards.mag, osample=wildcards.osample, fastq=fastqs)

def get_bai_osample(wildcards):
    fastqs = get_fastq_chunks(wildcards.osample)
    return expand(rules.minimap2_fq.output.bai, mag=wildcards.mag, osample=wildcards.osample, fastq=fastqs)

rule minimap2:
    input:
        bam = get_bam_osample,
        bai = get_bai_osample
    output:
        bam = os.path.join(results_path, "minimap2/{osample}/{mag}.bam"),
        bai = os.path.join(results_path, "minimap2/{osample}/{mag}.bam.bai"),
        stat = os.path.join(results_path, "minimap2/{osample}/{mag}.stat"),
    threads: 8
    log: os.path.join(results_path, "log/minimap2/{osample}/{mag}.log")
    shell:
        """
        module load minimap2 samtools

        lst=$(echo {output.bam} | sed "s/bam/list/")
        echo > $lst
        for bam in {input.bam} ; do
            echo $bam >> $lst
        done

        samtools merge -@ {threads} -o {output.bam} -b $lst --verbosity 1 2> {log}

        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.stat}
        """

def get_bam_osample_indv(wildcards):
    return rules.minimap2.output.bam

def get_bai_osample_indv(wildcards):
    return rules.minimap2.output.bai

rule coverm_per_sample:
    input:
        bam = get_bam_osample_indv,
        bai = get_bai_osample_indv
    output:
        os.path.join(results_path, "coverm/indv/{osample}/{mag}.tsv")
    conda: "coverm_"
    threads: 4
    log: os.path.join(results_path, "log/coverm_metacoag/{osample}/{mag}.log")
    message: "Coverm mini for {wildcards.osample}"
    shell:
        """
        export TMPDIR=$(dirname {output})"_tmp"
        mkdir -p $TMPDIR

        coverm contig -b {input.bam} \
         -p minimap2-ont \
         -o {output} \
         --min-read-percent-identity 95 \
         --contig-end-exclusion 0 \
         -t {threads} 2> {log}

        rm -rf $TMPDIR
        """

rule cm:
    input:
        bam = get_bam_osample_indv,
        bai = get_bai_osample_indv
    output:
        tmp = os.path.join(results_path, "coverm/metabat_method/{osample}/{mag}.tmp"),
        tsv = os.path.join(results_path, "coverm/metabat_method/{osample}/{mag}.tsv")
    conda: "coverm_"
    threads: 4
    log: os.path.join(results_path, "log/cm/{osample}/{mag}.log")
    message: "Coverm mini metabat for {wildcards.osample} mag {wildcards.mag}"
    params:
        pident = config["pident"]
    shell:
        """
        export TMPDIR=$(dirname {output.tsv})"_tmp"
        mkdir -p $TMPDIR

        coverm contig -v -b {input.bam} \
         -p minimap2-ont \
         -m metabat \
         -o {output.tmp} \
         --min-read-percent-identity {params.pident} \
         --contig-end-exclusion 0 \
         -t {threads} 2> /dev/null

        rm -rf $TMPDIR

        python src/generate_summary_per_mag_table.py -i {output.tmp} -o {output.tsv} -s {wildcards.osample} -m {wildcards.mag}
        """

def input_merge_metabat2_tables(wildcards):
    return expand(rules.cm.output.tmp, osample=get_samples_for_mag(wildcards), mag=wildcards.mag)
    # return expand(rules.cm.output.tmp, osample=samples_list, mag=wildcards.mag)

rule merge_metabat2_tables:
    input:
        input_merge_metabat2_tables
    output:
        real = os.path.join(results_path, "coverm/merge_metabat2_tables/{mag}.tsv"),
        check = os.path.join(results_path, "coverm/merge_metabat2_tables/{mag}.check")
    log:
        os.path.join(results_path, "log/merge_metabat2_tables/{mag}.log")
    priority: 10
    shell:
        """
        python src/merge_metabat2_tables_v2.py -i {input} -o {output.real} 2> {log}
        python src/merge_metabat2_tables_v2.py -i {input} -o {output.check} --keep-outliers
        """

rule magpurify2_coverage:
    input:
        coverage = rules.merge_metabat2_tables.output.check,
        genome = config["dereplicated_genome_path"]
    output:
        transform = os.path.join(results_path, "coverm/merge_metabat2_tables/{mag}.transform"),
        real = os.path.join(results_path, "magpurify2_coverage/{mag}/scores/coverage_scores.tsv")
    log:
        os.path.join(results_path, "log/magpurify2_coverage/{mag}.log")
    conda: "magpurify2_"
    params:
        min_identity = 0.95
    retries: 3
    shell:
        """
        outdir=$(dirname $(dirname {output.real}))
        rm -rf $outdir

        cat {input.coverage} | \
         ./src/transpose_table.py | \
         grep -v "var" | \
         grep -v "contigLen" | \
         grep -v "totalAvgDepth" | \
         ./src/transpose_table.py > {output.transform}

        magpurify2 coverage --coverage_file {output.transform} \
         --min_identity {params.min_identity} \
         {input.genome} $outdir 2> {log}
        """
        
def input_cm_all_samples_per_mag(wildcards):
    return expand(rules.cm.output.tsv, osample=samples_list, mag=wildcards.mag)

rule cm_all_samples_per_mag:
    input:
        input_cm_all_samples_per_mag
    output:
        os.path.join(results_path, "coverm/per_mag/{mag}.tsv"),
    shell:
        """
        python src/merge_2x2_tables.py -t {input} -o {output}
        """

rule cm_all:
    input:
        expand(rules.cm_all_samples_per_mag.output, mag=mags)
    output:
        os.path.join(results_path, "coverm/cm_all.tsv")
    shell:
        """
        python src/concat_tables.py -i {input} -o {output}
        """
