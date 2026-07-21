# QC report generation
# --------------------

rule fastqc:
    input:
        "data/trimmed/{sample}/{sample}_{bc}.fastq.gz"
    output:
       html = "data/fastqc/{sample}_{bc}_fastqc.html",
        zip = "data/fastqc/{sample}_{bc}_fastqc.zip"
    message:
        """--- Checking fastq files with FastQC."""
    params:
        outdir=lambda wc: f"data/fastqc/tmp_{wc.sample}_{wc.bc}",
        prefix=lambda wc: f"{wc.sample}_{wc.bc}"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc {input} -o {params.outdir}
        mv {params.outdir}/{params.prefix}_fastqc.html {output.html}
        mv {params.outdir}/{params.prefix}_fastqc.zip {output.zip}
        rm -rf {params.outdir}
        """


# run multiQC on output logs
# -----------------------------------------------------
rule multiqc:
    input:
        [f"data/fastqc/{sample}_{bc}_fastqc.zip"
        for sample in SAMPLES
        for bc in get_keep_list(sample)]
    output:
        report="results/multiqc_report.html"
    message:
        """--- Generating MultiQC report."""
    params:
        outdir="results",
        name="multiqc_report.html"
    shell:
        "multiqc -n {params.name} -o {params.outdir} {input}"