#############################################
# edgeR full analysis + QC plots (pre- and post-normalization)
# Author: Shashanka Puranika K
# Date: 24-11-2025
# Changes:
#  - fixed blocky smoothScatter by increasing grid and PNG resolution
#  - removed use of limma::decideTestsD on glmQLFTest object; MD plot now uses FDR-based status
#  - small safety checks to avoid 'plot.new has not been called yet' errors
#############################################

# --------- Config ---------
counts_file <- "gene_counts_clean.matrix.tsv"   # your input counts file
design_file <- "design.txt"                     # design file with 'Sample' and 'Group' columns
outdir <- "plots_edger"                         # where plots & html go
results_csv <- "edgeR_results.csv"              # DE results output
top_n_heatmap <- 50                             # top DE genes in heatmap

# Normalization options:
# - norm_method: "TMM" (default edgeR), or "RLE", or "upperquartile", or "none"
# - use_voom: if TRUE, use limma::voom() transformed values for downstream plots (v$E)
norm_method <- "TMM"
use_voom <- FALSE

# Volcano plot parameters
volcano_pCutoff <- 0.05
volcano_FCcutoff <- 1
volcano_top_n_labels <- 10   # number of genes to label on the static volcano
fdr_cutoff <- 0.05          # FDR threshold used for highlighting MD/Smear plots

# --------- Helpers ---------
install_if_missing <- function(pkgs) {
  to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(to_install)) {
    message("Installing missing packages: ", paste(to_install, collapse = ", "))
    install.packages(to_install, repos = "https://cloud.r-project.org")
  }
}

# Install/load packages
install_if_missing(c("edgeR", "pheatmap", "limma", "ggplot2", "ggrepel", "viridis", "dplyr"))
suppressPackageStartupMessages({
  library(edgeR)
  library(pheatmap)
  library(limma)     # used for voom if requested and some plotting helpers
  library(ggplot2)
  library(ggrepel)
  library(viridis)
  library(dplyr)
})

# --------- I/O checks ---------
if (!file.exists(counts_file)) {
  stop(sprintf("Counts file not found: %s\nSet 'counts_file' to the correct path.", counts_file))
}
if (!file.exists(design_file)) {
  stop(sprintf("Design file not found: %s\nSet 'design_file' to the correct path.", design_file))
}
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# --------- Data load ---------
counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
cat("Counts matrix loaded:\n")
print(head(counts))

# Read design file (must have Sample, Group)
design_df <- read.table(design_file, header = TRUE, sep = "", stringsAsFactors = FALSE, comment.char = "")
stopifnot(all(c("Sample","Group") %in% colnames(design_df)))
design_df$Sample <- trimws(design_df$Sample)
design_df$Group  <- trimws(design_df$Group)

# Align counts to design order
count_samples <- colnames(counts)
if (!all(design_df$Sample %in% count_samples)) {
  missing <- setdiff(design_df$Sample, count_samples)
  stop(sprintf("Counts matrix is missing %d samples listed in design.txt: %s",
               length(missing), paste(missing, collapse=", ")))
}
extra <- setdiff(count_samples, design_df$Sample)
if (length(extra)) {
  warning("Counts has ", length(extra), " samples not present in design.txt; they will be dropped: ",
          paste(extra, collapse=", "))
}
counts <- counts[, design_df$Sample, drop = FALSE]

# Build group factor & design matrix
group <- factor(tolower(design_df$Group))
levels(group) <- make.names(levels(group))
design <- model.matrix(~ 0 + group)
colnames(design) <- sub("^group", "", colnames(design))
cat("\nDesign matrix (first rows):\n")
print(design[1:min(6, nrow(design)), , drop=FALSE])

# --------- edgeR setup & filtering (keeps genes) ---------
grp_levels <- colnames(design)
if (length(grp_levels) == 2) {
  contrast_str <- sprintf("%s - %s", grp_levels[2], grp_levels[1])
  message("Using contrast: ", contrast_str)
} else {
  stop("Design has ", length(grp_levels), " groups. Please modify the script to define your desired contrast.")
}

dge <- DGEList(counts = counts, group = group)

# Filter lowly-expressed genes
keep <- filterByExpr(dge, design = design)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Save a pre-normalization snapshot (norm.factors = 1 by default)
pre_dge <- dge
# ensure pre_dge has norm.factors = 1
pre_dge$samples$norm.factors <- rep(1, nrow(pre_dge$samples))

# Pre-normalization logCPM (uses raw library sizes, no norm factors)
pre_logCPM <- cpm(pre_dge, log = TRUE, prior.count = 2, normalized.lib.sizes = FALSE)

# --------- PLOTS BEFORE NORMALIZATION ---------
# File naming: prefix "pre_"
# 1) Library sizes (raw)
png(file.path(outdir, "pre_01_library_sizes_barplot.png"), width = 1200, height = 800)
barplot(pre_dge$samples$lib.size / 1e6, names = colnames(pre_dge),
        las = 2, col = "steelblue", ylab = "Library size (millions)",
        main = "Pre-normalization: library sizes per sample")
abline(h = mean(pre_dge$samples$lib.size / 1e6), lty = 2, col = "red")
dev.off()

# 2) Density of raw counts (post-filter)
png(file.path(outdir, "pre_02_log10_counts_density.png"), width = 1200, height = 800)
plot(density(log10(pre_dge$counts[,1] + 1)), col="steelblue", lwd=2,
     main="Pre-normalization: Density of log10 raw counts (post-filter)", xlab="log10(count + 1)")
if (ncol(pre_dge$counts) > 1) {
  for (i in 2:ncol(pre_dge$counts)) lines(density(log10(pre_dge$counts[,i] + 1)), col=i+1, lwd=2)
}
legend("topright", legend=colnames(pre_dge), col=2:(ncol(pre_dge$counts)+1), lwd=2, bty="n")
dev.off()

# 3) Boxplot of pre-normalization logCPM
png(file.path(outdir, "pre_03_logCPM_boxplot.png"), width = 1200, height = 800)
boxplot(as.data.frame(pre_logCPM), las=2, col="gray80",
        main="Pre-normalization: Boxplot of logCPM (no library normalization)", ylab="logCPM")
dev.off()

# 4) MDS
png(file.path(outdir, "pre_04_MDS_plot.png"), width = 1200, height = 1000)
plotMDS(pre_dge, labels = colnames(pre_dge), col = as.numeric(pre_dge$samples$group),
        main="Pre-normalization: MDS (leading logFC) - samples")
legend("topright", legend=levels(pre_dge$samples$group),
       col=1:length(levels(pre_dge$samples$group)), pch=16, bty="n")
dev.off()

# 5) Sample correlation heatmap (pre)
cors_pre <- cor(pre_logCPM, method = "pearson")
png(file.path(outdir, "pre_05_sample_correlation_heatmap.png"), width = 1200, height = 1200)
op <- par(mar=c(6,6,2,2))
image(1:ncol(cors_pre), 1:ncol(cors_pre), t(cors_pre[nrow(cors_pre):1, ]),
      axes=FALSE, col=colorRampPalette(c("navy","white","firebrick"))(100),
      main="Pre-normalization: Sample-sample Pearson correlation (logCPM)")
axis(1, at=1:ncol(cors_pre), labels=colnames(cors_pre), las=2)
axis(2, at=1:ncol(cors_pre), labels=rev(colnames(cors_pre)), las=2)
par(op)
dev.off()

# 6) PCA (pre)
pca_pre <- prcomp(t(pre_logCPM), scale. = TRUE)
pct_pre <- round(100 * pca_pre$sdev^2 / sum(pca_pre$sdev^2), 1)
png(file.path(outdir, "pre_06_PCA.png"), width = 1200, height = 900)
plot(pca_pre$x[,1], pca_pre$x[,2], pch=19, col=as.numeric(pre_dge$samples$group),
     xlab=paste0("PC1 (", pct_pre[1], "%)"), ylab=paste0("PC2 (", pct_pre[2], "%)"),
     main="Pre-normalization: PCA of samples (logCPM)")
text(pca_pre$x[,1], pca_pre$x[,2], labels=colnames(pre_dge), pos=3, cex=0.9)
legend("topright", legend=levels(pre_dge$samples$group), col=1:length(levels(pre_dge$samples$group)), pch=19, bty="n")
dev.off()

# 7) Mean–variance trend (pre)
meanCPM_pre <- rowMeans(cpm(pre_dge, log=FALSE, normalized.lib.sizes = FALSE))
varCPM_pre  <- apply(cpm(pre_dge, log=FALSE, normalized.lib.sizes = FALSE), 1, var)

# --- Updated smoothScatter: increase nbin and PNG resolution to avoid blockiness ---
png(file.path(outdir, "pre_07_mean_variance_trend_CPM.png"), width = 2400, height = 1800, res = 150)
smoothScatter(log10(meanCPM_pre + 1), log10(varCPM_pre + 1),
              nbin = c(512, 512),
              colramp = colorRampPalette(c("white", "lightblue", "skyblue", "dodgerblue", "navy")),
              nrpoints = 0,
              xlab="log10(mean CPM + 1)", ylab="log10(variance CPM + 1)",
              main="Pre-normalization: Mean–variance trend (CPM)")
set.seed(1)
samp_idx <- sample.int(length(meanCPM_pre), size = min(5000, length(meanCPM_pre)))
points(log10(meanCPM_pre[samp_idx] + 1), log10(varCPM_pre[samp_idx] + 1),
       pch = 20, cex = 0.45, col = rgb(0,0,0,0.25))
dev.off()

# --------- NORMALIZATION ---------
if (!norm_method %in% c("TMM", "RLE", "upperquartile", "none")) {
  stop("Invalid norm_method. Choose one of 'TMM', 'RLE', 'upperquartile', or 'none'.")
}

if (norm_method == "none") {
  dge$samples$norm.factors <- rep(1, nrow(dge$samples))
  message("Normalization skipped (norm_method = 'none').")
} else {
  dge <- calcNormFactors(dge, method = norm_method)
  message("Applied normalization method: ", norm_method)
}

# Write sample table with normalization factors
sample_info <- dge$samples
write.csv(sample_info, file = file.path(outdir, "samples_with_libsizes_and_normfactors.csv"), row.names = TRUE)

# post-normalization logCPM (normalized.lib.sizes = TRUE)
post_logCPM <- if (use_voom) {
  voom_res <- voom(dge, design = design, plot = FALSE)
  voom_res$E
} else {
  cpm(dge, log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE)
}
# Also save normalized CPM (linear)
normCPM <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
write.csv(normCPM, file = file.path(outdir, "normalized_CPM.csv"), row.names = TRUE)
if (use_voom) write.csv(post_logCPM, file = file.path(outdir, "voom_E_logCPM.csv"), row.names = TRUE)

# --------- PLOTS AFTER NORMALIZATION ---------
# File naming: prefix "post_"
# 1) Effective library sizes barplot (post)
effective_lib_sizes_post <- dge$samples$lib.size * dge$samples$norm.factors
png(file.path(outdir, "post_01_library_sizes_barplot.png"), width = 1200, height = 800)
barplot(effective_lib_sizes_post/1e6, names = colnames(dge),
        las = 2, col = "steelblue", ylab = "Effective library size (millions)",
        main = paste0("Post-normalization: effective library sizes per sample (", norm_method, ")"))
abline(h = mean(effective_lib_sizes_post/1e6), lty = 2, col = "red")
dev.off()

# 2) Density of counts (post) - we still plot raw counts density, normalization affects effective libs
png(file.path(outdir, "post_02_log10_counts_density.png"), width = 1200, height = 800)
plot(density(log10(dge$counts[,1] + 1)), col="steelblue", lwd=2,
     main="Post-normalization: Density of log10 raw counts (post-filter)", xlab="log10(count + 1)")
if (ncol(dge$counts) > 1) {
  for (i in 2:ncol(dge$counts)) lines(density(log10(dge$counts[,i] + 1)), col=i+1, lwd=2)
}
legend("topright", legend=colnames(dge), col=2:(ncol(dge$counts)+1), lwd=2, bty="n")
dev.off()

# 3) Boxplot of post-normalization logCPM
png(file.path(outdir, "post_03_logCPM_boxplot.png"), width = 1200, height = 800)
boxplot(as.data.frame(post_logCPM), las=2, col="gray80",
        main=paste0("Post-normalization: Boxplot of logCPM (", norm_method, ifelse(use_voom, " + voom", ""), ")"),
        ylab="logCPM")
dev.off()

# 4) MDS (post)
png(file.path(outdir, "post_04_MDS_plot.png"), width = 1200, height = 1000)
plotMDS(dge, labels = colnames(dge), col = as.numeric(dge$samples$group),
        main="Post-normalization: MDS (leading logFC) - samples")
legend("topright", legend=levels(dge$samples$group),
       col=1:length(levels(dge$samples$group)), pch=16, bty="n")
dev.off()

# 5) Sample correlation heatmap (post)
cors_post <- cor(post_logCPM, method = "pearson")
png(file.path(outdir, "post_05_sample_correlation_heatmap.png"), width = 1200, height = 1200)
op <- par(mar=c(6,6,2,2))
image(1:ncol(cors_post), 1:ncol(cors_post), t(cors_post[nrow(cors_post):1, ]),
      axes=FALSE, col=colorRampPalette(c("navy","white","firebrick"))(100),
      main="Post-normalization: Sample-sample Pearson correlation (logCPM)")
axis(1, at=1:ncol(cors_post), labels=colnames(cors_post), las=2)
axis(2, at=1:ncol(cors_post), labels=rev(colnames(cors_post)), las=2)
par(op)
dev.off()

# 6) PCA (post)
pca_post <- prcomp(t(post_logCPM), scale. = TRUE)
pct_post <- round(100 * pca_post$sdev^2 / sum(pca_post$sdev^2), 1)
png(file.path(outdir, "post_06_PCA.png"), width = 1200, height = 900)
plot(pca_post$x[,1], pca_post$x[,2], pch=19, col=as.numeric(dge$samples$group),
     xlab=paste0("PC1 (", pct_post[1], "%)"), ylab=paste0("PC2 (", pct_post[2], "%)"),
     main="Post-normalization: PCA of samples (logCPM)")
text(pca_post$x[,1], pca_post$x[,2], labels=colnames(dge), pos=3, cex=0.9)
legend("topright", legend=levels(dge$samples$group), col=1:length(levels(dge$samples$group)), pch=19, bty="n")
dev.off()

# 7) Mean–variance trend (post)
meanCPM_post <- rowMeans(cpm(dge, log=FALSE, normalized.lib.sizes = TRUE))
varCPM_post  <- apply(cpm(dge, log=FALSE, normalized.lib.sizes = TRUE), 1, var)

# --- Updated smoothScatter for post-normalization as well ---
png(file.path(outdir, "post_07_mean_variance_trend_CPM.png"), width = 2400, height = 1800, res = 150)
smoothScatter(log10(meanCPM_post + 1), log10(varCPM_post + 1),
              nbin = c(512, 512),
              colramp = colorRampPalette(c("white", "lightblue", "skyblue", "dodgerblue", "navy")),
              nrpoints = 0,
              xlab="log10(mean CPM + 1)", ylab="log10(variance CPM + 1)",
              main="Post-normalization: Mean–variance trend (normalized CPM)")
set.seed(1)
samp_idx2 <- sample.int(length(meanCPM_post), size = min(5000, length(meanCPM_post)))
points(log10(meanCPM_post[samp_idx2] + 1), log10(varCPM_post[samp_idx2] + 1),
       pch = 20, cex = 0.45, col = rgb(0,0,0,0.25))
dev.off()

# --------- edgeR differential testing (after normalization) ---------
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
contrast <- makeContrasts(contrasts = contrast_str, levels = design)
qlf <- glmQLFTest(fit, contrast = contrast)

# Show top tags
cat("\nTopTags (head):\n")
print(topTags(qlf))

# Full results
results <- topTags(qlf, n = nrow(dge$counts))$table
write.csv(results, file = results_csv, row.names = TRUE)
cat(sprintf("\nedgeR results written to: %s\n", results_csv))

# --------- DE-based visualizations (post-normalization) ---------
# Static volcano with ggplot2 + ggrepel
results_df <- as.data.frame(results)
results_df$Gene <- rownames(results_df)
res_plot <- results_df %>%
  mutate(negLogP = -log10(pmax(PValue, .Machine$double.xmin)),
         absLogFC = abs(logFC),
         signif = "NS") %>%
  mutate(signif = case_when(
    absLogFC >= volcano_FCcutoff & PValue < volcano_pCutoff ~ "P & Log2FC",
    PValue < volcano_pCutoff ~ "P",
    absLogFC >= volcano_FCcutoff ~ "Log2FC",
    TRUE ~ "NS"
  ))

volcano_colors <- c("NS" = "grey30", "Log2FC" = "royalblue", "P" = "forestgreen", "P & Log2FC" = "red3")
label_genes <- res_plot %>% arrange(FDR, desc(absLogFC)) %>% head(volcano_top_n_labels) %>% pull(Gene)

p_volcano_post <- ggplot(res_plot, aes(x = logFC, y = negLogP, color = signif)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(values = volcano_colors) +
  geom_vline(xintercept = c(-volcano_FCcutoff, volcano_FCcutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(volcano_pCutoff), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = paste("Post-normalization Volcano:", contrast_str),
       subtitle = paste0("P cutoff = ", volcano_pCutoff, "; |log2FC| cutoff = ", volcano_FCcutoff),
       x = "log2 Fold Change", y = "-log10(P-value)", color = "") +
  theme(legend.position = "right") +
  geom_text_repel(data = res_plot %>% filter(Gene %in% label_genes),
                  aes(label = Gene), size = 3, max.overlaps = 20)

png(file.path(outdir, "post_08_volcano_ggplot.png"), width = 1600, height = 1200)
print(p_volcano_post)
dev.off()

# MD plot (post)
# Instead of decideTestsD (which may not take glmQLFTest objects), build a logical status vector
# using the FDR from the topTags results and pass it to plotMD.
md_status <- NULL
if (!is.null(results$FDR)) {
  md_status <- results$FDR < fdr_cutoff
  # ensure the status vector is named and matches the qlf table order if possible
  names(md_status) <- rownames(results)
}

png(file.path(outdir, "post_09_MD_plot.png"), width = 1200, height = 900)
# Use tryCatch so that abline is only called if plotting succeeded.
tryCatch({
  # plotMD will highlight points where md_status is TRUE
  if (!is.null(md_status)) {
    # pass the logical status (named) to highlight significant genes
    plotMD(qlf, status = md_status, column = 1, cex = 0.6,
           main = "Post-normalization: MD plot (contrast)", xlab = "Average logCPM", ylab = "LogFC")
  } else {
    plotMD(qlf, column = 1, cex = 0.6,
           main = "Post-normalization: MD plot (contrast)", xlab = "Average logCPM", ylab = "LogFC")
  }
  abline(h = 0, col = "gray40", lty = 2)
}, error = function(e) {
  message("MD plot failed: ", e$message)
})
dev.off()

# MA plot (post)
avg_logCPM_post <- rowMeans(post_logCPM)
png(file.path(outdir, "post_10_MA_plot.png"), width = 1200, height = 900)
plot(avg_logCPM_post, results$logFC, pch=20, col=rgb(0,0,0,0.4),
     xlab="Average logCPM", ylab="Log2 Fold Change", main="Post-normalization: MA plot")
abline(h=0, col="red", lty=2)
dev.off()

# plotSmear (post)
png(file.path(outdir, "post_11_smear_plot.png"), width = 1200, height = 900)
# highlight de.tags using FDR threshold
de_tags <- rownames(results)[results$FDR < fdr_cutoff]
if (length(de_tags) == 0) de_tags <- NULL
plotSmear(qlf, de.tags = de_tags)
abline(h = c(-1, 1), col = "blue")
dev.off()

# p-value histogram (post)
png(file.path(outdir, "post_12_pvalue_histogram.png"), width = 1200, height = 900)
hist(results$PValue, breaks=50, col="gray70", border="white",
     main="Post-normalization: Histogram of raw p-values", xlab="P-value")
dev.off()

# Top genes heatmap (post)
top_genes <- rownames(results[order(results$FDR), ])[1:min(top_n_heatmap, nrow(results))]
mat <- post_logCPM[top_genes, , drop = FALSE]
mat_scaled <- t(scale(t(mat)))
ann <- data.frame(Group = dge$samples$group)
rownames(ann) <- colnames(dge)

png(file.path(outdir, "post_13_topDE_heatmap.png"), width = 1600, height = 1600)
pheatmap(mat_scaled,
         annotation_col = ann,
         show_rownames = TRUE, show_colnames = TRUE,
         color = colorRampPalette(c("navy","white","firebrick"))(100),
         main = paste0("Post-normalization: Top ", min(top_n_heatmap, nrow(results)), " DE genes (z-scored logCPM)"))
dev.off()

# --------- HTML index & session info ---------
html_file <- file.path(outdir, "index.html")
sink(html_file)
cat("<html><head><title>edgeR QC & DE plots (pre/post normalization)</title></head><body>\n")
cat("<h1>edgeR QC & DE plots (pre- and post-normalization)</h1>\n")
imgs <- list.files(outdir, pattern="\\.png$", full.names = FALSE)
for (f in imgs) {
  cat(sprintf("<div><h3>%s</h3><img src=\"%s\" style=\"max-width:100%%;height:auto;\"></div>\n", f, f))
}
cat("</body></html>\n")
sink()

writeLines(capture.output(sessionInfo()),
           con = file.path(outdir, "sessionInfo.txt"))

cat(sprintf("\nAll plots saved in: %s\nBrowse via: %s\n", outdir, html_file))
#############################################
# End of script
#############################################
