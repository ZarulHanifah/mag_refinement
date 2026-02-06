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

rule longstitch:
    input:
        mag = config["dereplicated_genome_path"],
        fq = lambda wc: os.path.join(temp_path, f"minimap2_fq_collected/{wc.mag}/all.fq")
    output:
        # tmp_dir = temp(directory(os.path.join(results_path, "longstitch/parking/{mag}"))),
        # tmp_mag = temp(os.path.join(results_path, "longstitch/parking/{mag}/{mag}.fa")),
        mag = os.path.join(results_path, "longstitch/{mag}/{mag}.longstitch.fasta"),
        supp = directory(os.path.join(results_path, "longstitch/{mag}/other_files"))
    log:
        os.path.join(results_path, "log/longstitch/{mag}.log")
    conda: "longstitch_"
    params:
        gsize = 5_000_000
    threads: 4
    shell:
        """
        outdir=$(dirname {output.mag})
        out_prefix=$outdir"/"{wildcards.mag}
        tmp_dir=$(dirname $outdir)"/parking/"{wildcards.mag}

        tmp_mag=$tmp_dir"/"draft".fa"
        tmp_fq=$tmp_dir"/"reads".fq"

        mkdir -p $tmp_dir

        ./src/mlFASTA2slFASTA.sh < $(readlink -f {input.mag} ) > $tmp_mag
        ln -sf {input.fq} $tmp_fq

        cd $tmp_dir

        longstitch run draft=draft \
         reads=reads \
         G={params.gsize} \
         t={threads}  &> {log}

        cp $(find $tmp_dir -type f | grep abyss-scaffold.fa) {output.mag}
        mkdir -p {output.supp}
        cp $(find $tmp_dir -type f | grep tigmint | grep -v "fa$") {output.supp}

        cd -
        rm -rf $tmp_dir
        """

