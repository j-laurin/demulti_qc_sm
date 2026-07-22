# Demultiplexing of raw reads
# -------------------------

rule demultiplex:
    input:
        fastq="data/{sample}.fastq.gz",
        barcodes="config/barcodes.fa"
    output:
        temp(expand("data/trimmed/{{sample}}/{sample}_{barcode}.fastq.gz", barcode=BARCODES, sample=SAMPLES)),
        done=touch("logs/{sample}_3_umi_cut.done")
    message:
        """--- Demultiplexing based on adaptor sequences"""
    threads: 8
    params:
        outdir = lambda wc: f"data/trimmed/{wc.sample}",
        length = config["length_cutoff"]
    log:
        step1 = "logs/{sample}_1_adapter.log",
        step2 = "logs/{sample}_2_demux.log"

    # 1 - cutting constant adaptor with max 2 mismatches, no indels, 
    # keeping reads longer than 10 nt and trimmed
    # 2 - demultiplexing by the barcodes.fa file in config folder, 
    # keeping reads with barcode/index and longer than 10 nt
    # 3 - removing 3 nt from 5'(start) and 6 nt from 3'(end), 
    # keeping reads longer than 19 nt

    shell:
        """
        mkdir -p {params.outdir}/tmp_demux logs/{wildcards.sample}_3_umi_cut

        cutadapt -j {threads} --no-indels -q 20,20 --trimmed-only \
            -a CTGTAGGCACCATCAAT -e 0.15 -m 10 {input.fastq} 2> {log.step1} | \
        cutadapt -j {threads} --no-indels --trimmed-only \
            -a file:{input.barcodes} -m 10 -o {params.outdir}/tmp_demux/{{name}}.fastq.gz - > {log.step2} 2>&1

        for f in {params.outdir}/tmp_demux/*.fastq.gz; do
            bc=$(basename "$f" .fastq.gz)
            cutadapt -j {threads} -u 3 -u -6 -m {params.length} \
            -o {params.outdir}/{wildcards.sample}_${{bc}}.fastq.gz "$f" \
            > logs/{wildcards.sample}_3_umi_cut/${{bc}}.log 2>&1
        done

        rm -rf {params.outdir}/tmp_demux
        """


N_LIBRARIES = 8

rule cutadapt_stats:
    input:
        step1 = expand("logs/{sample}_1_adapter.log", sample=SAMPLES),
        demux = expand("logs/{sample}_2_demux.log", sample=SAMPLES),
        # step3_dirs = expand("logs/{sample}_3_umi_cut", sample=SAMPLES),
        flag = expand("logs/{sample}_3_umi_cut.done", sample=SAMPLES),
    output:
        pdf = "results/read_distribution_plots.pdf",
        csv = "results/read_distribution_table.csv"
    params:
        samples=SAMPLES,
        n_libraries=N_LIBRARIES,
        step3_dirs = expand("logs/{sample}_3_umi_cut", sample=SAMPLES)
    script:
        "../scripts/reads_progress_plots.py"