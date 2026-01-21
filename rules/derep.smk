

def get_dereplicated_genomes_path(wildcards):
    return config["dereplicated_genomes"][wildcards.sample]

rule copy_genomes:
    input:
        get_dereplicated_genomes_path
    output:
        done = touch(os.path.join(results_path_ie79, "tmp/copy_genomes/{sample}.done"))
    shell:
        """
        outdir=$(dirname {output.done})
        mkdir -p $outdir

        for i in $(ls {input} | grep "fasta$") ; do
            j2=$outdir"/{wildcards.sample}_"$i
            cp {input}/$i $j2
        done
        """

rule drep:
    input:
        expand(rules.copy_genomes.output, sample=samples_list)
    output:
        clust = os.path.join(results_path_ie79, "derep/data_tables/Cdb.csv"),
        paths = os.path.join(results_path_ie79, "derep/data_tables/Bdb.csv"),
        data  = directory(os.path.join(results_path_ie79, "derep/data"))
    conda:   "drep_"
    threads: 12
    log:     os.path.join(results_path_ie79, "log/drep/log.log")
    message: "Clustering MAGs"
    shell:
        """
        indir=$(dirname {input[0]})
        outdir=$(dirname $(dirname {output.clust}))
        rm -rf $outdir

        dRep compare $outdir -p {threads} -g $indir/*fasta 2> {log}
        """

summary_list = [ config["summary"][sample] for sample in samples_list ]

rule concat_summary:
    input:
        summary_list
    output:
        os.path.join(results_path_ie79, "concat_summary.tsv")
    shell:
        """
        head -1 {input[0]} > {output}

        for i in {input}; do
            sample=$(basename $(dirname $i))
            sample=$(echo $sample | sed "s/00//" | sed "s/_//")
            cat $i | sed "1d" | sed "s/^/$sample\_/" >> {output}
        done
        """

rule compile_drep:
    input:
        clust  = rules.drep.output.clust,
        paths  = rules.drep.output.paths,
        checkm = rules.concat_summary.output
    output:
        outdir = directory(os.path.join(results_path_ie79, "derep/dereplicated_genomes")),
        summa  = os.path.join(results_path_ie79, "derep/checkm2_summary.tsv"),
        summa2 = "summary.tsv"
    log: os.path.join(results_path_ie79, "log/compile_drep/log.log")
    shell:
        """
        python src/process_cluster_checkm.py --cluster {input.clust} \
         --checkm {input.checkm} \
         -o {output.summa} -l {log}

        mkdir -p {output.outdir}

        for i in $(cat {output.summa} | cut -f1) ; do
            p=$(grep -w $i".fasta" {input.paths} | tr "," "\t" | cut -f2)
            cp -r $p {output.outdir}
        done

        head -1 {input.checkm} > {output.summa2}

        for i in $(ls {output.outdir} ; do
         grep -w $i {input.checkm}
        done | sort -k 12,12 >> {output.summa2}
        """
