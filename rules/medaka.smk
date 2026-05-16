
rule index_original_genome:
    input: config["dereplicated_genome_path"],
    output: os.path.join(results_path, "resources/index_original_genome/{mag}/assem.mmi")
    conda: "medaka_"
    log: os.path.join(results_path, "log/index_original_genome/{mag}.log")
    threads: 8
    shell:
        """
        minimap2 -x map-ont -d {output} {input} 2> {log}
        """

rule minialign:
    input:
        fq = lambda wc: os.path.join(temp_path, f"minimap2_fq_collected/{wc.mag}/all.fq"),
        assem = config["dereplicated_genome_path"],
        assem_index = rules.index_original_genome.output
    output:
        bam = temp(os.path.join(results_path, "resources/minialign/{mag}.bam")),
        bai = temp(os.path.join(results_path, "resources/minialign/{mag}.bam.bai")),
    conda: "medaka_"
    log:
        os.path.join(results_path, "log/minialign/{mag}.log")
    threads: 2
    shell:
        """
        prefix=$(echo {output.bam} | sed "s/\.bam$//")

        mini_align -i {input.fq} \
         -r {input.assem} \
         -m -p $prefix \
         -t {threads} 2> {log}

        samtools index {output.bam}
        """

rule medaka_inference:
    input:
        bam = rules.minialign.output.bam,
        bai = rules.minialign.output.bai,
    output:
        os.path.join(results_path, "medaka_inference/{mag}.hdf")
    conda: "medaka_"
    log:
        os.path.join(results_path, "log/medaka_inference/{mag}.log")
    threads: 2
    params:
        model = medaka_model
    shell:
        """
        medaka inference \
         --model {params.model} \
         --threads {threads} \
         {input.bam} {output} &> {log}
        """

rule medaka:
    input:
        hdf = rules.medaka_inference.output,
        assem = config["dereplicated_genome_path"]
    output:
        assem = os.path.join(results_path, "medaka/{mag}.fasta")
    log:
        os.path.join(results_path, "log/medaka/{mag}.log")
    conda: "medaka_"
    threads: 2
    shell:
        """
        medaka sequence {input.hdf} {input.assem} {output} 2> {log}
        rm -rf {output}".gaps_in_draft_coords.bed"
        """

rule proovframe_o6:
    input:
        assem = rules.medaka.output.assem,
        uniref_db = uniref_db
    output:
        temp(os.path.join(results_path, "proovframe/o6/{mag}.o6"))
    conda: "proovframe_"
    log: os.path.join(results_path, "log/proovframe/{mag}.log")
    threads: 20
    shell:
        """
        proovframe map -t {threads} \
         -d {input.uniref_db} \
         -o {output} \
         -m fast \
         {input.assem} 2> {log} 
        """

rule proovframe:
    input:
        assem = rules.medaka.output.assem,
        o6 = rules.proovframe_o6.output
    output:
        assem = os.path.join(results_path, "proovframe/assem/{mag}.fasta")
    conda: "proovframe_"
    log: os.path.join(results_path, "log/proovframe/{mag}.log")
    shell:
        """
        proovframe fix -D -o {output} \
         {input.assem} {input.o6} 2> {log}
        """
# medaka round 2

rule index_medaka_rd2:
    input: rules.medaka.output.assem
    output: os.path.join(results_path, "resources/index_medaka_rd2/{mag}/assem.mmi")
    conda: "medaka_"
    log: os.path.join(results_path, "log/index_medaka_rd2/{mag}.log")
    threads: 8
    shell:
        """
        minimap2 -x map-ont -d {output} {input} 2> {log}
        """

rule minialign_rd2:
    input:
        fq = lambda wc: os.path.join(temp_path, f"minimap2_fq_collected/{wc.mag}/all.fq"),
        assem = rules.medaka.output.assem,
        assem_index = rules.index_medaka_rd2.output
    output:
        bam = temp(os.path.join(results_path, "resources/minialign_rd2/{mag}.bam")),
        bai = temp(os.path.join(results_path, "resources/minialign_rd2/{mag}.bam.bai")),
    conda: "medaka_"
    log: os.path.join(results_path, "log/minialign_rd2/{mag}.log")
    threads: 2
    shell:
        """
        prefix=$(echo {output.bam} | sed "s/\.bam$//")

        mini_align -i {input.fq} \
         -r {input.assem} \
         -m -p $prefix \
         -t {threads} 2> {log}

        samtools index {output.bam}
        """

rule medaka_inference_rd2:
    input:
        bam = rules.minialign_rd2.output.bam,
        bai = rules.minialign_rd2.output.bai,
    output:
        os.path.join(results_path, "medaka_inference_rd2/{mag}.hdf")
    conda: "medaka_"
    log:
        os.path.join(results_path, "log/medaka_inference_rd2/{mag}.log")
    threads: 2
    params:
        model = medaka_model
    shell:
        """
        medaka inference \
         --model {params.model} \
         --threads {threads} \
         {input.bam} {output} &> {log}
        """

rule medaka_rd2:
    input:
        hdf = rules.medaka_inference_rd2.output,
        assem = rules.medaka.output.assem
    output:
        assem = os.path.join(results_path, "medaka_rd2/{mag}.fasta")
    log:
        os.path.join(results_path, "log/medaka_rd2/{mag}.log")
    conda: "medaka_"
    threads: 2
    shell:
        """
        medaka sequence {input.hdf} {input.assem} {output} 2> {log}
        rm -rf {output}".gaps_in_draft_coords.bed"
        """
