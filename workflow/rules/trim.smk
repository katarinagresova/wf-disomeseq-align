# Read preprocessing: adapter trim -> collapse identical reads (the random nt
# are UMIs) -> remove the random nt.


rule cutadapt_reads:
    input:
        lambda wc: samples.loc[wc.sample, "fastq_1"],
    output:
        temp(f"{RESULTS_DIR}/cutadapt_reads/{{sample}}.fastq.gz"),
    # cutadapt's report; rule multiqc and qc_summary read it
    log:
        f"{LOG_DIR}/cutadapt/{{sample}}.log",
    params:
        adapter=config["adapter"],
        # an insert of only random nt (adapter dimer) is removed here; the
        # footprint length filter is in trim_reads
        minlength=config["randomer_5p"] + config["randomer_3p"] + 1,
    conda:
        "../envs/cutadapt.yaml"
    threads: 4
    shell:
        # -O 6: default -O 3 lets a 3 nt chance match to the adapter trim real
        # insert bases off the 3' end; -O 6 is long enough to be a real adapter.
        # --discard-untrimmed: a footprint is shorter than the read, so a read
        # without the adapter is not a footprint.
        "cutadapt --cores {threads} -a {params.adapter} -O 6 --discard-untrimmed "
        "--minimum-length {params.minlength} -o {output} {input} > {log} 2>&1"


# Collapse BEFORE removing the random nt: they are part of the sequence here,
# so identical footprints with different random nt stay distinct molecules.
# The log has the read counts and lengths (with the random nt) for qc_summary.
rule collapse_reads:
    input:
        fastq=f"{RESULTS_DIR}/cutadapt_reads/{{sample}}.fastq.gz",
        script=workflow.source_path("../scripts/collapse_reads.py"),
    output:
        temp(f"{RESULTS_DIR}/collapse_reads/{{sample}}.fastq.gz"),
    log:
        f"{LOG_DIR}/collapse_reads/{{sample}}.log",
    conda:
        "../envs/pysam.yaml"
    # holds every distinct read in memory
    resources:
        mem_mb=16000,
    shell:
        "python {input.script} {input.fastq} {wildcards.sample} 2> {log} | gzip > {output}"


rule trim_reads:
    input:
        f"{RESULTS_DIR}/collapse_reads/{{sample}}.fastq.gz",
    output:
        f"{RESULTS_DIR}/trim_reads/{{sample}}.fastq.gz",
    log:
        f"{LOG_DIR}/trim_reads/{{sample}}.log",
    params:
        randomer_5p=config["randomer_5p"],
        randomer_3p=config["randomer_3p"],
        minlength=config["minlength"],
    conda:
        "../envs/cutadapt.yaml"
    threads: 4
    shell:
        # Cut the random nt off both ends of each read and keep them in the read
        # name, e.g. (4 + 4) read "s_0_x12" that starts with ACCG and ends with
        # TATA becomes "s_0_x12_ACCG:TATA" (x12: 12 identical reads).
        "cutadapt --cores {threads} -u {params.randomer_5p} -u -{params.randomer_3p} "
        "--minimum-length {params.minlength} "
        "--rename '{{header}}_{{cut_prefix}}:{{cut_suffix}}' -o {output} {input} > {log} 2>&1"
