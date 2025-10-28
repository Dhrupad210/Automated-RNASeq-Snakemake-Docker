library(DESeq2)
library(ggplot2)
library(ggrepel)

# Get config
cfg <- yaml::yaml.load_file("config.yaml")

# Extract conditions
sample_names <- names(cfg$samples)
conditions <- as.factor(unname(cfg$samples))

# Read count tables
count_files <- paste0("results/", sample_names, ".counts.txt")
count_list <- lapply(count_files, function(f) read.table(f, header=TRUE, row.names=1, comment.char="#")[, "Count"])
counts <- do.call(cbind, count_list)
colnames(counts) <- sample_names

# DESeq2 analysis
coldata <- data.frame(row.names=sample_names, condition=conditions)
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)
dds <- DESeq(dds)
res <- results(dds)

# Save results
write.csv(as.data.frame(res), "results/deseq2_results.csv")

# Volcano plot
res$gene <- rownames(res)
res$significant <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "Significant", "Not Significant")

ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=significant, label=ifelse(padj < 0.05 & abs(log2FoldChange) > 1, gene, ""))) +
  geom_point() +
  geom_text_repel(size=3) +
  theme_minimal() +
  labs(title="Volcano Plot", x="log2 Fold Change", y="-log10 Adjusted p-value") +
  scale_color_manual(values=c("grey", "red")) +
  theme(plot.title=element_text(hjust=0.5))

ggsave("results/volcano_plot.png")

