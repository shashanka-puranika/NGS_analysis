#############################################
# edgeR full analysis + QC plots using design.txt
# Author: Shashanka Puranika K
# Date: 24-11-2025
# Modifications: replaced plotly interactive volcano with ggiraph-based interactive volcano
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
volcano_top_n_labels <- 10   # number of genes to label on the static volcano and interactive

# --------- Helpers ---------
install_if_missing <- function(pkgs) {
  to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(to_install)) {
    message("Installing missing packages: ", paste(to_install, collapse = ", "))
    install.packages(to_install, repos = "https://cloud.r-project.org")
  }
}

# Install/load packages
install_if_missing(c("edgeR", "pheatmap", "limma", "ggplot2", "ggrepel", "ggiraph", "htmlwidgets", "viridis", "dplyr"))
suppressPackageStartupMessages({
  library(edgeR)
  library(pheatmap)
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(ggiraph)
  library(htmlwidgets)
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

# Read design file (must have Sample, Group). Example:
# Sample  Group
# SRR4420293 control
# SRR4420294 control
# SRR4420297 test
# SRR4420298 test
design_df <- read.table(design_file, header = TRUE, sep = "", stringsAsFactors = FALSE, comment.char = "")
stopifnot(all(c("Sample","Group") %in% colnames(design_df)))  # requires columns named 'Sample' and 'Group'
design_df$Sample <- trimws(design_df$Sample)
design_df$Group  <- trimws(design_df$Group)

# Align counts to design order
count_samples <- colnames(counts)
if (!all(design_df$Sample %in% count_samples)) {
  missing <- setdiff(design_df$Sample, count_samples)
  stop(sprintf("Counts matrix is missing %d samples listed in design.txt: %s",
               length(missing), paste(missing, collapse=", ")))
}
# Warn if counts has extra columns not in design
extra <- setdiff(count_samples, design_df$Sample)
if (length(extra)) {
  warning("Counts has ", length(extra), " samples not present in design.txt; they will be dropped: ",
          paste(extra, collapse=", "))
}
# Reorder / subset columns to design order
counts <- counts[, design_df$Sample, drop = FALSE]

# Build group factor from design (lowercase & safe names)
group <- factor(tolower(design_df$Group))
levels(group) <- make.names(levels(group))  # safe column names for model.matrix

# Design matrix with no intercept
design <- model.matrix(~ 0 + group)
colnames(design) <- sub("^group", "", colnames(design))
cat("\nDesign matrix (first rows):\n")
print(design[1:min(6, nrow(design)), , drop=FALSE])

# --------- edgeR pipeline with automatic contrast ---------
grp_levels <- colnames(design)
if (length(grp_levels) == 2) {
  contrast_str <- sprintf("%s - %s", grp_levels[2], grp_levels[1])
  message("Using contrast: ", contrast_str)  # e.g., 'test - control'
} else {
  stop("Design has ", length(grp_levels), " groups. Please modify the script to define your desired contrast (e.g., case - control).")
}

dge <- DGEList(counts = counts, group = group)

# Keep genes with reasonable counts
keep <- filterByExpr(dge, design = design)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# --------- NORMALIZATION ---------
# Acceptable methods for calcNormFactors: "TMM", "RLE", "upperquartile"
if (!norm_method %in% c("TMM", "RLE", "upperquartile", "none")) {
  stop("Invalid norm_method. Choose one of 'TMM', 'RLE', 'upperquartile', or 'none'.")
}

if (norm_method == "none") {
  # set all normalization factors to 1 (no normalization)
  dge$samples$norm.factors <- rep(1, nrow(dge$samples))
  message("Normalization skipped (norm_method = 'none').")
} else {
  dge <- calcNormFactors(dge, method = norm_method)
  message("Applied normalization method: ", norm_method)
}

# Write sample table with normalization factors
sample_info <- dge$samples
write.csv(sample_info, file = file.path(outdir, "samples_with_libsizes_and_normfactors.csv"), row.names = TRUE)

# Normalized counts as CPM (counts per million, normalized library sizes)
normCPM <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
write.csv(normCPM, file = file.path(outdir, "normalized_CPM.csv"), row.names = TRUE)
# Also produce logCPM used for plots (either voom or edgeR log-CPM)
if (use_voom) {
  voom_res <- voom(dge, design = design, plot = FALSE)
  logCPM <- voom_res$E       # log2 CPM-ish values from voom
  write.csv(logCPM, file = file.path(outdir, "voom_E_logCPM.csv"), row.names = TRUE)
  message("Using voom-transformed values for downstream plots (use_voom = TRUE).")
} else {
  logCPM <- cpm(dge, log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE)
  message("Using edgeR logCPM (prior.count = 2) from normalized CPM for downstream plots.")
}

# --------- Continue edgeR differential testing ---------
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

contrast <- makeContrasts(contrasts = contrast_str, levels = design)
qlf <- glmQLFTest(fit, contrast = contrast)

# Show top tags in console
cat("\nTopTags (head):\n")
print(topTags(qlf))

# Full results
results <- topTags(qlf, n = nrow(dge$counts))$table
write.csv(results, file = results_csv, row.names = TRUE)
cat(sprintf("\nedgeR results written to: %s\n", results_csv))

# --------- Volcano using ggplot2 + ggrepel + ggiraph (replaces plotly) ---------
# Prepare results for plotting
results_df <- as.data.frame(results)
results_df$Gene <- rownames(results_df)

res_plot <- results_df %>%
  dplyr::mutate(negLogP = -log10(pmax(PValue, .Machine$double.xmin)),
                absLogFC = abs(logFC),
                signif = "NS") %>%
  dplyr::mutate(signif = dplyr::case_when(
    absLogFC >= volcano_FCcutoff & PValue < volcano_pCutoff ~ "P & Log2FC",
    PValue < volcano_pCutoff ~ "P",
    absLogFC >= volcano_FCcutoff ~ "Log2FC",
    TRUE ~ "NS"
  ))

# Color mapping
volcano_colors <- c("NS" = "grey30", "Log2FC" = "royalblue", "P" = "forestgreen", "P & Log2FC" = "red3")

# Which genes to label: order by FDR then by abs(logFC)
label_genes <- res_plot %>% arrange(FDR, desc(absLogFC)) %>% head(volcano_top_n_labels) %>% pull(Gene)

# Static ggplot2 volcano (PNG)
p_volcano <- ggplot(res_plot, aes(x = logFC, y = negLogP, color = signif)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(values = volcano_colors) +
  geom_vline(xintercept = c(-volcano_FCcutoff, volcano_FCcutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(volcano_pCutoff), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = paste("Volcano plot:", contrast_str),
       subtitle = paste0("P cutoff = ", volcano_pCutoff, "; |log2FC| cutoff = ", volcano_FCcutoff),
       x = "log2 Fold Change", y = "-log10(P-value)", color = "") +
  theme(legend.position = "right") +
  geom_text_repel(data = res_plot %>% filter(Gene %in% label_genes),
                  aes(label = Gene),
                  size = 3, max.overlaps = 20)

png(file.path(outdir, "00_volcano_ggplot.png"), width = 1600, height = 1200)
print(p_volcano)
dev.off()

# Interactive volcano via ggiraph, saved as standalone HTML
# Add tooltip text
res_plot <- res_plot %>%
  dplyr::mutate(tooltip = paste0(Gene, "\nlogFC=", signif(logFC, 3),
                                 "\nP=", signif(PValue, 3),
                                 "\nFDR=", signif(FDR, 3)))

# interactive ggplot: use geom_point_interactive and geom_text_interactive for labeled genes
p_volcano_int <- ggplot(res_plot, aes(x = logFC, y = negLogP, color = signif)) +
  ggiraph::geom_point_interactive(aes(tooltip = tooltip, data_id = Gene), alpha = 0.7, size = 2) +
  scale_color_manual(values = volcano_colors) +
  geom_vline(xintercept = c(-volcano_FCcutoff, volcano_FCcutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(volcano_pCutoff), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = paste("Interactive Volcano:", contrast_str),
       subtitle = paste0("P cutoff = ", volcano_pCutoff, "; |log2FC| cutoff = ", volcano_FCcutoff),
       x = "log2 Fold Change", y = "-log10(P-value)", color = "") +
  theme(legend.position = "right")

# add interactive labels for the top genes
label_df <- res_plot %>% filter(Gene %in% label_genes)
if (nrow(label_df) > 0) {
  p_volcano_int <- p_volcano_int +
    ggiraph::geom_text_interactive(data = label_df,
                                   aes(label = Gene, tooltip = tooltip, data_id = Gene),
                                   size = 3,
                                   vjust = -1)
}

girafe_obj <- girafe(ggobj = p_volcano_int,
                     options = list(
                       opts_tooltip(css = "background-color:white;padding:5px;border:1px solid #333;border-radius:4px;"),
                       opts_hover(css = "stroke-width:2;")
                     ))

htmlwidgets::saveWidget(girafe_obj, file = file.path(outdir, "00_volcano_interactive.html"), selfcontained = TRUE)

# --------- QC & Exploratory Plots ---------

# 1) Basic QC: library sizes & distributions
effective_lib_sizes <- dge$samples$lib.size * dge$samples$norm.factors
png(file.path(outdir, "01_library_sizes_barplot.png"), width = 1200, height = 800)
barplot(effective_lib_sizes/1e6, names = colnames(dge),
        las = 2, col = "steelblue", ylab = "Effective library size (millions)",
        main = paste0("Effective library sizes per sample (method: ", norm_method, ")"))
abline(h = mean(effective_lib_sizes/1e6), lty = 2, col = "red")
dev.off()

png(file.path(outdir, "02_log10_counts_density.png"), width = 1200, height = 800)
plot(density(log10(dge$counts[,1] + 1)), col="steelblue", lwd=2,
     main="Density of log10 raw counts (post-filter)", xlab="log10(count + 1)")
if (ncol(dge$counts) > 1) {
  for (i in 2:ncol(dge$counts)) lines(density(log10(dge$counts[,i] + 1)), col=i+1, lwd=2)
}
legend("topright", legend=colnames(dge), col=2:(ncol(dge$counts)+1), lwd=2, bty="n")
dev.off()

png(file.path(outdir, "03_logCPM_boxplot.png"), width = 1200, height = 800)
boxplot(as.data.frame(logCPM), las=2, col="gray80",
        main=paste0("Boxplot of logCPM (normalized, method: ", norm_method, ifelse(use_voom, " + voom", ""), ")"),
        ylab="logCPM")
dev.off()

# 2) Sample relationships: MDS, correlation heatmap, PCA
png(file.path(outdir, "04_MDS_plot.png"), width = 1200, height = 1000)
plotMDS(dge, labels = colnames(dge), col = as.numeric(dge$samples$group),
        main="MDS (leading logFC) - samples")
legend("topright", legend=levels(dge$samples$group),
       col=1:length(levels(dge$samples$group)), pch=16, bty="n")
dev.off()

cors <- cor(logCPM, method = "pearson")
png(file.path(outdir, "05_sample_correlation_heatmap.png"), width = 1200, height = 1200)
op <- par(mar=c(6,6,2,2))
image(1:ncol(cors), 1:ncol(cors), t(cors[nrow(cors):1, ]),
      axes=FALSE, col=colorRampPalette(c("navy","white","firebrick"))(100),
      main="Sample-sample Pearson correlation (logCPM)")
axis(1, at=1:ncol(cors), labels=colnames(cors), las=2)
axis(2, at=1:ncol(cors), labels=rev(colnames(cors)), las=2)
par(op)
dev.off()

pca <- prcomp(t(logCPM), scale. = TRUE)
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
groups <- dge$samples$group
cols <- as.numeric(groups)
png(file.path(outdir, "06_PCA.png"), width = 1200, height = 900)
plot(pca$x[,1], pca$x[,2], pch=19, col=cols,
     xlab=paste0("PC1 (", pct[1], "%)"), ylab=paste0("PC2 (", pct[2], "%)"),
     main="PCA of samples (logCPM)")
text(pca$x[,1], pca$x[,2], labels=colnames(dge), pos=3, cex=0.9)
legend("topright", legend=levels(groups), col=1:length(levels(groups)), pch=19, bty="n")
dev.off()

# 3) Dispersion diagnostics
png(file.path(outdir, "07_BCV_plot.png"), width = 1200, height = 900)
plotBCV(dge, main="BCV plot after estimateDisp")
dev.off()

png(file.path(outdir, "08_QL_dispersion_plot.png"), width = 1200, height = 900)
plotQLDisp(fit, main="Quasi-likelihood dispersion")
dev.off()

# 4) DE signal & effect-size visualizations
png(file.path(outdir, "09_MD_plot.png"), width = 1200, height = 900)
plotMD(qlf, status = decideTestsD(qlf), column=1, cex=0.6,
       main="MD plot: Contrast", xlab="Average logCPM", ylab="LogFC")
abline(h = 0, col = "gray40", lty = 2)
dev.off()

avg_logCPM <- rowMeans(logCPM)
png(file.path(outdir, "10_MA_plot.png"), width = 1200, height = 900)
plot(avg_logCPM, results$logFC, pch=20, col=rgb(0,0,0,0.4),
     xlab="Average logCPM", ylab="Log2 Fold Change", main="MA plot")
abline(h=0, col="red", lty=2)
dev.off()

png(file.path(outdir, "11_smear_plot.png"), width = 1200, height = 900)
plotSmear(qlf, de.tags = rownames(results)[results$FDR < 0.05])
abline(h = c(-1, 1), col = "blue")
dev.off()

png(file.path(outdir, "12_pvalue_histogram.png"), width = 1200, height = 900)
hist(results$PValue, breaks=50, col="gray70", border="white",
     main="Histogram of raw p-values", xlab="P-value")
dev.off()

# 5) Top genes heatmap
top_genes <- rownames(results[order(results$FDR), ])[1:min(top_n_heatmap, nrow(results))]
mat <- logCPM[top_genes, ]
mat_scaled <- t(scale(t(mat)))
ann <- data.frame(Group = dge$samples$group)
rownames(ann) <- colnames(dge)

png(file.path(outdir, "14_topDE_heatmap.png"), width = 1600, height = 1600)
pheatmap(mat_scaled,
         annotation_col = ann,
         show_rownames = TRUE, show_colnames = TRUE,
         color = colorRampPalette(c("navy","white","firebrick"))(100),
         main = paste0("Top ", min(top_n_heatmap, nrow(results)), " DE genes (z-scored logCPM)"))
dev.off()

# 6) Mean–variance trend
meanCPM <- rowMeans(cpm(dge, log=FALSE, normalized.lib.sizes = TRUE))
varCPM  <- apply(cpm(dge, log=FALSE, normalized.lib.sizes = TRUE), 1, var)

png(file.path(outdir, "15_mean_variance_trend_CPM.png"), width = 1200, height = 900)
smoothScatter(log10(meanCPM + 1), log10(varCPM + 1),
              xlab="log10(mean CPM + 1)", ylab="log10(variance CPM + 1)",
              main="Mean–variance trend (normalized CPM)")
dev.off()

# --------- HTML index & session info ---------
html_file <- file.path(outdir, "index.html")
sink(html_file)
cat("<html><head><title>edgeR QC & DE plots</title></head><body>\n")
cat("<h1>edgeR QC & DE plots</h1>\n")
# include PNGs and HTML interactive plots
files_to_show <- list.files(outdir, pattern="\\.(png|html)$", full.names = TRUE)
for (f in files_to_show) {
  fname <- basename(f)
  if (grepl("\\.png$", f)) {
    cat(sprintf("<div><h3>%s</h3><img src=\"%s\" style=\"max-width:100%%;height:auto;\"></div>\n", fname, fname))
  } else {
    cat(sprintf("<div><h3>%s</h3><a href=\"%s\">Open interactive plot</a></div>\n", fname, fname))
  }
}
cat("</body></html>\n")
sink()

writeLines(capture.output(sessionInfo()),
           con = file.path(outdir, "sessionInfo.txt"))

cat(sprintf("\nAll plots saved in: %s\nBrowse via: %s\n", outdir, html_file))
#############################################
# End of script
#############################################
