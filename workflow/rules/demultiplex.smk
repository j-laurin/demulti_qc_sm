# Demultiplexing of raw reads
# -------------------------

rule demultiplex:
    input:
        fastq="data/{sample}.fastq.gz",
        barcodes="config/barcodes.fa"
    output:
        temp(expand("data/trimmed/{{sample}}/{barcode}.fastq.gz", barcode=BARCODES))
    log:
        "logs/{sample}.txt"
    message:
        """--- Demultiplexing based on adaptor sequences"""
    threads: 8
    params:
        outdir=lambda wc: f"data/trimmed/{wc.sample}"
    shell:
        """
        mkdir -p {params.outdir}
        cutadapt -g file:{input.barcodes} \
            -o {params.outdir}/{{name}}.fastq.gz \
            {input.fastq} > {params.outdir}/cutadapt.log
        """
