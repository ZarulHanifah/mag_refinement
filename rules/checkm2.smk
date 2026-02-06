rule checkm2_original:
    input:
        mag = config["dereplicated_genome_path"],
        db = checkm2_db
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

