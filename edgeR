###############################################################################
# Differential Gene Expression Analysis using edgeR
# Input  : featureCounts gene_counts.txt
# Output : DEG tables, plots, normalized counts
#
# Experimental design example:
#   4 samples total
#   2 Control vs 2 Treated
#
# Author : Standard Academic RNA-seq Pipeline
###############################################################################

############################
# STEP 0: Load libraries
############################
suppressMessages({
  library(edgeR)
  library(ggplot2)
  library(pheatmap)
  library(data.table)
})

############################
# STEP 1: Read featureCounts output
############################
# featureCounts output has:
#  - Comment lines starting with #
#  - GeneID in column 1
#  - Counts start from column 7 onward

counts_raw <- fread(
  "counts/gene_counts.txt",
  data.table = FALSE
)

# Remove annotation columns, keep only count matrix
count_matrix <- as.matrix(counts_raw[, 7:ncol(counts_raw)])
rownames(count_matrix) <- counts_raw$Geneid

# Check dimensions
cat("Genes:", nrow(count_matrix), "\n")
cat("Samples:", ncol(count_matrix), "\n")

############################
# STEP 2: Define sample metadata
############################
# IMPORTANT: Order must match columns in count_matrix

sample_info <- data.frame(
  sample = colnames(count_matrix),
  condition = factor(c(
    "Control", "Control",
    "Treated", "Treated"
  ))
)

print(sample_info)

############################
# STEP 3: Create edgeR DGEList object
############################
dge <- DGEList(
  counts = count_matrix,
  group  = sample_info$condition
)

############################
# STEP 4: Filter lowly expressed genes
############################
# edgeR recommendation:
# Keep genes expressed in enough samples

keep_genes <- filterByExpr(dge)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

cat("Genes after filtering:", nrow(dge), "\n")

############################
# STEP 5: Library size normalization
############################
# Uses TMM (Trimmed Mean of M-values)

dge <- calcNormFactors(dge)

# View normalization factors
dge$samples

############################
# STEP 6: Experimental design matrix
############################
design <- model.matrix(~ condition, data = sample_info)
colnames(design)

############################
# STEP 7: Estimate dispersion
############################
# Biological variability estimation
# Required before differential testing

dge <- estimateDisp(dge, design)

# Visualize dispersion
plotBCV(dge)

############################
# STEP 8: Fit GLM model
############################
fit <- glmQLFit(dge, design)

############################
# STEP 9: Differential expression testing
############################
# Coefficient 2 = conditionTreated vs Control

qlf <- glmQLFTest(fit, coef = 2)

# Extract results table
deg_table <- topTags(qlf, n = Inf)$table

# Add FDR explicitly (already included, but clarity matters)
deg_table$FDR <- p.adjust(deg_table$PValue, method = "BH")

############################
# STEP 10: Save DEG results
############################
write.csv(
  deg_table,
  file = "results/edgeR_DEG_full_results.csv",
  row.names = TRUE
)

############################
# STEP 11: Define significant DEGs
############################
sig_deg <- deg_table[
  abs(deg_table$logFC) >= 1 &
  deg_table$FDR < 0.05,
]

cat("Significant DEGs:", nrow(sig_deg), "\n")

write.csv(
  sig_deg,
  file = "results/edgeR_significant_DEGs.csv",
  row.names = TRUE
)

############################
# STEP 12: Volcano plot
############################
volcano_data <- data.frame(
  logFC = deg_table$logFC,
  negLogFDR = -log10(deg_table$FDR)
)

ggplot(volcano_data, aes(logFC, negLogFDR)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Volcano Plot (edgeR)",
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  )

############################
# STEP 13: MDS plot (sample similarity)
############################
plotMDS(dge,
        labels = sample_info$condition,
        col = c("blue", "blue", "red", "red"))

############################
# STEP 14: Heatmap of top DEGs
############################
top_genes <- rownames(sig_deg)[1:min(50, nrow(sig_deg))]

logCPM <- cpm(dge, log = TRUE)

pheatmap(
  logCPM[top_genes, ],
  scale = "row",
  annotation_col = sample_info["condition"],
  show_rownames = TRUE,
  fontsize_row = 8,
  main = "Top Differentially Expressed Genes"
)

############################
# STEP 15: MA plot
############################
plotMD(qlf)
abline(h = c(-1, 1), col = "blue")

###############################################################################
# END OF edgeR PIPELINE
###############################################################################
