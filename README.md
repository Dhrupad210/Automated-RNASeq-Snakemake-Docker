# Automated-RNASeq-Snakemake-Docker
A reproducible RNA-Seq analysis pipeline built with Snakemake, R, and Docker, designed for easy setup and cross-platform reproducibility.

# 🧬 RNA-Seq Analysis Pipeline (Snakemake + Docker)

A reproducible RNA-Seq analysis pipeline built with Snakemake, R, and Docker, designed for easy setup and cross-platform reproducibility.

The pipeline performs:

1. Quality control and trimming (fastp)

2. Alignment (hisat2)

3. BAM processing (samtools)

4. Quantification (featureCounts)

5. Differential expression analysis (DESeq2)

6. Volcano plot visualization with gene labels (ggplot2 + ggrepel)

📁 Folder Structure
rna_seq_pipeline/

├── Snakefile · config.yaml · Dockerfile
├── scripts/ → deseq2_analysis.R
├── data/ → control1.fastq.gz, control2.fastq.gz, treated1.fastq.gz, treated2.fastq.gz
├── reference/ → genome.fa, annotation.gtf

└── results/



🐍 Snakefile

Save this as Snakefile:

```
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
  ```

🧾 config.yaml

Save as config.yaml:
```
samples:
  control1: control
  control2: control
  treated1: treated
  treated2: treated
```

📜 deseq2_analysis.R

Save as scripts/deseq2_analysis.R:
```
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(yaml)

cfg <- yaml::yaml.load_file("config.yaml")

sample_names <- names(cfg$samples)
conditions <- as.factor(unname(cfg$samples))

count_files <- paste0("results/", sample_names, ".counts.txt")
count_list <- lapply(count_files, function(f) read.table(f, header=TRUE, row.names=1, comment.char="#")[, "Count"])
counts <- do.call(cbind, count_list)
colnames(counts) <- sample_names

coldata <- data.frame(row.names=sample_names, condition=conditions)
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
dds <- DESeq(dds)
res <- results(dds)

write.csv(as.data.frame(res), "results/deseq2_results.csv")

res$gene <- rownames(res)
res$significant <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "Significant", "Not Significant")

ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=significant,
                label=ifelse(padj < 0.05 & abs(log2FoldChange) > 1, gene, ""))) +
  geom_point() +
  geom_text_repel(size=3) +
  theme_minimal() +
  labs(title="Volcano Plot", x="log2 Fold Change", y="-log10 Adjusted p-value") +
  scale_color_manual(values=c("grey", "red")) +
  theme(plot.title=element_text(hjust=0.5))

ggsave("results/volcano_plot.png")
```

🐳 Dockerfile

Save this as Dockerfile:
```
FROM continuumio/miniconda3:latest

LABEL maintainer="Dhrupad Banerjee <245hsbb010@ibab.ac.in>"
LABEL description="Dockerized RNA-Seq Snakemake Pipeline"

WORKDIR /pipeline

RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    graphviz \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y -c bioconda -c conda-forge \
    snakemake \
    fastp \
    hisat2 \
    samtools \
    subread \
    bioconductor-deseq2 \
    r-ggplot2 \
    r-ggrepel \
    r-yaml

COPY . /pipeline

ENTRYPOINT ["snakemake"]
CMD ["--help"]
```

🧱 Build the Docker Image
`docker build -t rna-seq-pipeline`

▶️ Run the Pipeline
```
docker run -it --rm \
  -v $(pwd):/pipeline \
  -w /pipeline \
  rna-seq-pipeline \
  --cores 4
```

✅ Results will appear in the results/ folder on your machine.

🧪 Test Run (Dry Mode)

To check workflow structure without running commands:

```
docker run -it --rm \
  -v $(pwd):/pipeline \
  -w /pipeline \
  rna-seq-pipeline \
  --dry-run
```

📊 Output Files
results/
├── control1.counts.txt

├── control2.counts.txt

├── treated1.counts.txt

├── treated2.counts.txt

├── deseq2_results.csv

└── volcano_plot.png


The Volcano plot highlights significantly up/downregulated genes (adjusted p < 0.05 and |log₂FC| > 1).

🧠 Tips & Notes

To visualize rule dependencies:

` snakemake --dag | dot -Tpng > dag.png `


To rerun a single rule:

`snakemake results/deseq2_results.csv --cores 4`


Edit config.yaml to add or remove samples — Snakemake will adapt automatically.
