rule checkm2_original:
    input:
        mag = ancient(config["dereplicated_genome_path"]),
        db = ancient(checkm2_db)
    output:
        remove_list = temp([
            directory(os.path.join(results_path, "checkm2/original/{mag}/protein_files/"))
        ]),
        tmp = temp(os.path.join(results_path, ".tmp/checkm2/original/{mag}/{mag}.fasta")),
        report = os.path.join(results_path, "checkm2/original/{mag}/quality_report.tsv")
    conda: "checkm2_"
    threads: 2
    log: os.path.join(results_path, "log/checkm2/original/{mag}.log")
    message: "Running checkm2 for genome {wildcards.mag}: original"
    shell:
        """
        cp {input.mag} {output.tmp}
        input_genome=$(find {output.tmp} | grep "fasta")

        outfolder=$(dirname {output.report})

        checkm2 predict --force \
         --threads {threads} \
         --input $input_genome \
         --database_path {input.db} \
         --output-directory $outfolder &> {log}
        """

rule checkm1_original:
    input:
        mag = ancient(config["dereplicated_genome_path"]),
        db = ancient(checkm1_db)
    output:
        tmpdir = temp([ 
            directory(os.path.join(results_path, "checkm1_tmp/original/{mag}"))
        ]),
        outdir = directory(os.path.join(results_path, "checkm1/original/{mag}")) 
    log: os.path.join(results_path, "log/checkm1/original/{mag}.log")
    conda: "checkm_"
    threads: 2
    shell:
        """
        export CHECKM_DATA_PATH={input.db}
        module load hmmer

        mkdir -p {output.tmpdir} ; cp {input.mag} {output.tmpdir}

        outdir=$()

        checkm lineage_wf -t {threads} \
         -x fasta \
         {output.tmpdir} {output.outdir} &> {log}
        """
