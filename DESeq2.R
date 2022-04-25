library(DESeq2)

counts <- read.table("isomiR_counts.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
groups <- read.table("annotation.tsv", sep = "\t", header = TRUE)
expr_group <- "left"
ctrl_group <- "right"

groups$Tissue <- relevel(factor(groups$Tissue), ref = ctrl_group)
rownames(groups) <- groups$Sample
groups$Sample <- NULL

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = groups,
                              design = ~ Tissue)

dds <- DESeq(dds)
res <- lfcShrink(dds, coef = paste("Tissue_", expr_group, "_vs_", ctrl_group, sep=""), type = "apeglm")
res <- res[order(res$padj),]

write.table(res, file = "gene_deseq2_R.tsv", sep = "\t", quote = FALSE)
