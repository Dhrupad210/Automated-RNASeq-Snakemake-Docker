rule all:
    input:
        expand("results/{sample}.counts.txt", sample=["control1", "control2", "treated1", "treated2"]),
        "results/deseq2_results.csv",
        "results/volcano_plot.png"

rule fastp:
    input:
        "data/{sample}.fastq.gz"
    output:
        "results/{sample}_clean.fastq.gz"
    shell:
        "fastp -i {input} -o {output} -q 20 -l 50"

rule hisat2:
    input:
        reads = "results/{sample}_clean.fastq.gz",
        ref = "reference/genome.fa"
    output:
        "results/{sample}.sam"
    shell:
        "hisat2 -x {input.ref} -U {input.reads} -S {output}"

rule samtools_sort:
    input:
        "results/{sample}.sam"
    output:
        "results/{sample}.bam"
    shell:
        "samtools sort -o {output} {input}"

rule featurecounts:
    input:
        bam = "results/{sample}.bam",
        gtf = "reference/annotation.gtf"
    output:
        "results/{sample}.counts.txt"
    shell:
        "featureCounts -a {input.gtf} -o {output} {input.bam}"

rule deseq2:
    input:
        counts = expand("results/{sample}.counts.txt", sample=["control1", "control2", "treated1", "treated2"]),
        cfg = "config.yaml"
    output:
        "results/deseq2_results.csv",
        "results/volcano_plot.png"
    script:
        "scripts/deseq2_analysis.R"

