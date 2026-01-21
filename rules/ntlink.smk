rule ntlink_prelim:
    input:
        mag = config["dereplicated_genome_path"],
        fq = lambda wc: os.path.join(temp_path, f"minimap2_fq_collected/{wc.mag}/all.fq")
    output:
        directory(os.path.join(results_path, "ntlink_prelim/{mag}"))
    log:
        os.path.join(results_path, "log/ntlink_prelim/{mag}.log")
    conda: "ntlink_"
    params:
        rounds = 5
    threads: 4
    shell:
        """
        prefix={output}"/"{wildcards.mag}
        target={output}"/"{wildcards.mag}.fasta

        mkdir -p {output}
        cp {input.mag} {output}

        ntLink scaffold gap_fill target=$target reads={input.fq} \
         prefix=$prefix t={threads} &> {log}


        #ntLink_rounds run_rounds target={input.mag} reads={input.fq} \
        # rounds={params.rounds}
        """

