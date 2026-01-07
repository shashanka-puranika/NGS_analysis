library("limma")

targets <- readTargets("Targets.txt")

RG <- read.maimages(targets$FileName, source="agilent")

files <- targets[,c("Cy3","Cy5")]

x <- read.maimages(targets, source="agilent", green.only=TRUE)

summary(x$E)

plotMD(RG)

plotMD(RG, array=2, main="Sample 2")

plotMD(RG, array=3, main="Sample 3")

plotMD(RG, array=4, main="Sample 4")

plotMD(RG, array=5, main="Sample 5")

plotMD(RG, array=6, main="Sample 6")

boxplot(data.frame(log2(RG$Gb)),main="Green background")

boxplot(data.frame(log2(RG$Rb)),main="Red background")

imageplot(log2(RG$Gb[,1]),RG$printer)

plotDensities(RG)

MA.p <-normalizeWithinArrays(RG)

plotDensities(MA.p)

MA.pAq <- normalizeBetweenArrays(MA.p, method="Aquantile")

plotDensities(MA.pAq)

MA.q <- normalizeBetweenArrays(RG, method="quantile")

plotDensities(MA.q, col="blue")

design <- modelMatrix(targets, ref ="normal")
fit <- lmFit(MA.q,design)
eb <- eBayes(fit)
tb <- topTable (eb, adjust.method="BH")
tt <- topTable(eb, n=nrow(MA.q))


topGenes <- tt[tt[, "P.Value"]<0.05,]
select <- rownames(topGenes)
select <- as.numeric(select)
selectedProbes <- MA.q[select,]


exprs <- MA.q$M

gene_labels <- ifelse(is.na(MA.q$genes$GeneName) | MA.q$genes$GeneName == "",
                      
                      MA.q$genes$ProbeName,
                      
                      MA.q$genes$GeneName)
rownames(exprs) <- gene_labels
top50_idx <- head(order(tt$P.Value), 50)
expr_top50 <- exprs[top50_idx, ]
rownames(expr_top50) <- make.unique(rownames(expr_top50))
library(pheatmap)

pheatmap(expr_top50,
         
         scale = "row",
         
         clustering_distance_rows = "euclidean",
         
         clustering_distance_cols = "euclidean",
         
         clustering_method = "complete",
         
         show_rownames = TRUE,
         
         show_colnames = TRUE,
         
         main = "Top 50 Differentially Expressed Genes")       

library(tidyverse) 
library(cluster) 
library(factoextra) 
library(limma)
library(gridExtra)

# Differential Expression Analysis
targets <- readTargets("Targets.txt")
rg <- read.maimages(targets, source="agilent")
rg <- backgroundCorrect(rg, method="normexp")
ma <- normalizeWithinArrays(rg, method="loess")
ma <- normalizeBetweenArrays(ma, method="Aquantile")
design <- cbind(TumorVsNormal = rep(1, ncol(ma$M)))
colnames(design) <- "TumorVsNormal"
fit <- lmFit(ma, design)
fit <- eBayes(fit)
top_genes <- topTable(fit, n=50)

# K-means Clustering
df <- top_genes[, c("logFC", "P.Value")]
df <- scale(df)

# Determine optimal number of clusters
fviz_nbclust(df, kmeans, method = "wss")

# Perform K-means clustering for k = 2 to 7
k2 <- kmeans(df, centers = 2, nstart = 25)
k3 <- kmeans(df, centers = 3, nstart = 25)
k4 <- kmeans(df, centers = 4, nstart = 25)
k5 <- kmeans(df, centers = 5, nstart = 25)
k6 <- kmeans(df, centers = 6, nstart = 25)
k7 <- kmeans(df, centers = 7, nstart = 25)

# Visualize the clusters
p1 <- fviz_cluster(k2, geom = "point", data = df) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point", data = df) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point", data = df) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point", data = df) + ggtitle("k = 5")
p5 <- fviz_cluster(k6, geom = "point", data = df) + ggtitle("k = 6")
p6 <- fviz_cluster(k7, geom = "point", data = df) + ggtitle("k = 7")

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)

# Pairwise cluster plot for k = 2
df %>%
  as_tibble() %>%
  mutate(cluster = k2$cluster,
         gene = rownames(top_genes)) %>%
  ggplot(aes(logFC, P.Value, color = factor(cluster), label = gene)) +
  geom_text()
