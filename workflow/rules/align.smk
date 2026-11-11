# Align the trimmed reads to the transcripts (main species and spike-in) and
# the contaminants in one pass, and keep the reads with one best alignment, to
# a transcript: one BAM for each species.


rule bowtie_index:
    input:
        transcripts=expand(f"{RESULTS_DIR}/reference/{{species}}/transcripts.fa", species=SPECIES),
        contaminants=f"{RESULTS_DIR}/reference/contaminants.fa",
    output:
        multiext(f"{RESULTS_DIR}/bowtie_index/index", ".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt",
                 ".rev.2.ebwt"),
    log:
        f"{LOG_DIR}/bowtie_index.log",
    params:
        prefix=lambda wildcards, output: output[0].removesuffix(".1.ebwt"),
    conda:
        "../envs/bowtie.yaml"
    threads: 8
    resources:
        mem_mb=8000,
    shell:
        "bowtie-build --threads {threads} $(echo {input.transcripts} | tr ' ' ','),{input.contaminants} "
        "{params.prefix} > {log} 2>&1"


# BOWTIE_ARGS: forward strand only (R1 is the sense strand), end to end with at
# most 2 mismatches, up to 256 alignments per read, best first, in input order.
# assign_reads.py gives each read a class (see the script). Fails if bowtie
# skipped a read or the read counts differ.
BOWTIE_ARGS = "-v 2 --norc -k 256 --best --reorder"


rule bowtie_align:
    input:
        fastq=f"{RESULTS_DIR}/trim_reads/{{sample}}.fastq.gz",
        index=rules.bowtie_index.output,
        contaminants=f"{RESULTS_DIR}/reference/contaminants.fa",
        spike_in=expand(f"{RESULTS_DIR}/reference/{{species}}/transcripts.fa", species=SPECIES[1:]),
        script=workflow.source_path("../scripts/assign_reads.py"),
    output:
        bam=f"{RESULTS_DIR}/bowtie/main/{{sample}}.bam",
        bai=f"{RESULTS_DIR}/bowtie/main/{{sample}}.bam.bai",
        # no records without a spike-in
        spike_in_bam=f"{RESULTS_DIR}/bowtie/spike_in/{{sample}}.bam",
        spike_in_bai=f"{RESULTS_DIR}/bowtie/spike_in/{{sample}}.bam.bai",
        counts=f"{RESULTS_DIR}/bowtie/{{sample}}.counts.tsv",
    log:
        bowtie=f"{LOG_DIR}/bowtie/{{sample}}.log",
        assign=f"{LOG_DIR}/bowtie/{{sample}}.assign.log",
    params:
        index=lambda wildcards, input: input.index[0].removesuffix(".1.ebwt"),
        bowtie_args=BOWTIE_ARGS,
        spike_in=lambda wildcards, input: input.spike_in or "-",
    conda:
        "../envs/bowtie.yaml"
    threads: 8
    resources:
        mem_mb=8000,
    shell:
        r"""
        # the unsorted BAMs are not outputs: remove them also if a step fails
        trap 'rm -f {output.bam}.unsorted.bam {output.spike_in_bam}.unsorted.bam' EXIT
        bowtie -p {threads} {params.bowtie_args} --sam -x {params.index} {input.fastq} 2> {log.bowtie} |
            python {input.script} {input.contaminants} {params.spike_in} {output.bam}.unsorted.bam \
                {output.spike_in_bam}.unsorted.bam {output.counts} 2> {log.assign}
        if grep -q "Exhausted best-first chunk memory" {log.bowtie}; then
            echo "ERROR: bowtie skipped reads (Exhausted best-first chunk memory): raise --chunkmbs" >> {log.assign}
            exit 1
        fi
        processed=$(awk -F': ' '$1 == "# reads processed" {{print $2}}' {log.bowtie})
        assigned=$(awk -F'\t' '$1 == "molecules" {{print $2}}' {log.assign})
        if [ "$processed" != "$assigned" ]; then
            echo "ERROR: bowtie processed $processed reads, assign_reads.py assigned $assigned" >> {log.assign}
            exit 1
        fi
        for bam in {output.bam} {output.spike_in_bam}; do
            python -c 'import sys, pysam; pysam.sort("-@", sys.argv[1], "-m", "512M", "-o", sys.argv[2], sys.argv[3]); pysam.index(sys.argv[2])' \
                {threads} $bam $bam.unsorted.bam
            rm $bam.unsorted.bam
        done
        """
