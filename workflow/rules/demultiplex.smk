# Demultiplexing of raw reads
# -------------------------

rule demultiplex:
    input:
        fastq="data/{sample}.fastq.gz",
        barcodes="config/barcodes.fa"
    output:
        temp(expand("data/trimmed/{{sample}}/{sample}_{barcode}.fastq.gz", barcode=BARCODES, sample=SAMPLES))
    message:
        """--- Demultiplexing based on adaptor sequences"""
    threads: 8
    params:
        outdir=lambda wc: f"data/trimmed/{wc.sample}"
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
            cutadapt -j {threads} -u 3 -u -6 -m 19 \
            -o {params.outdir}/{wildcards.sample}_${{bc}}.fastq.gz "$f" \
            > logs/{wildcards.sample}_3_umi_cut/${{bc}}.log 2>&1
        done

        rm -rf {params.outdir}/tmp_demux
        """
