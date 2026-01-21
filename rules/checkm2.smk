rule checkm2_fq:
    input:
        fasta = rules.flye_fq.output.assem,
        db = config["checkm2_db"]
    output:
        tmp = temp(os.path.join(results_path, ".tmp/{mag}/{mag}.fasta")),
        report = os.path.join(results_path, "checkm2/flye_fq/{mag}/storage/quality_report.tsv")
    conda: "checkm2_"
    threads: 2
    log: os.path.join(results_path, "log/checkm2/flye_fq/{mag}.log")
    message: "Running checkm2 for genome {wildcards.mag}"
    shell:
        """
        cp {input.fasta} {output.tmp}
        input_genome=$(find {output.tmp} | grep "fasta")
        
        outfolder=$(dirname {output.report})

        checkm2 predict --threads {threads} \
         --input $input_genome \
         --database_path {input.db} \
         --output-directory $outfolder &> {log}
        """

# rule checkm2_herro:
#     input:
#         fasta = rules.flye_herro.output.assem,
#         db = config["checkm2_db"]
#     output:
#         tmp = temp(os.path.join(results_path, ".tmp/{mag}/{mag}.fasta")),
#         report = os.path.join(results_path, "checkm2/flye_herro/{mag}/storage/quality_report.tsv")
#     conda: "checkm2_"
#     threads: 2
#     log:
#         os.path.join(results_path, "log/checkm2/flye_herro/{mag}.log")
#     message: "Running checkm2 for genome {wildcards.mag}"
#     shell:
#         """
#         cp {input.fasta} {output.tmp}
#         input_genome=$(find {output.tmp} | grep "fasta")
#         
#         outfolder=$(dirname {output.report})
# 
#         checkm2 predict --threads {threads} \
#          --input $input_genome \
#          --database_path {input.db} \
#          --output-directory $outfolder &> {log}
#         """

