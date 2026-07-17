# Demultiplexing of raw reads
# -------------------------

rule demultiplex:
    input:
        fastq ="data/{sample}.fastq.gz",
        barcodes ="config/barcodes.fa"
    output:
        temp("results/{sample}.demulti.fastq.gz")
    log:
        "results/logs/{sample}.demulti.txt"
    message:
        """--- Demultiplexing based on adaptor sequences"""
    threads: 
        8
    shell:
        """
        (cutadapt -j {threads} -m 20 -O 20 -a "polyA=A{{20}}" -a "QUALITY=G{{20}}" -n 2 {input.fq} | \
        cutadapt -j {threads} -m 20 -O 3 --nextseq-trim=10 -a "r1adapter=A{{18}}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" - | \
        cutadapt -j {threads} -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" --discard-trimmed -o {output} -) \
        > {log} 2>&1
        """

