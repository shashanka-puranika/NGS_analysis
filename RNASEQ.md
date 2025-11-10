\# Create project directory structure

mkdir -p rnaseq\_analysis/{raw\_data,qc/{pre\_trim,post\_trim},trimmed,reference,alignment,counts,results}

cd rnaseq\_analysis



\# Install required tools (if not already installed)

\# FastQC, Trimmomatic, HISAT2, SAMtools, featureCounts, DESeq2

\# Additionally: gffread for GFF to GTF conversion



cd raw\_data



\# Control samples

wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/003/SRR4420293/SRR4420293\_1.fastq.gz

wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/003/SRR4420293/SRR4420293\_2.fastq.gz

wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/004/SRR4420294/SRR4420294\_1.fastq.gz

wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/004/SRR4420294/SRR4420294\_2.fastq.gz



\# Test/Treatment samples

wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/007/SRR4420297/SRR4420297\_1.fastq.gz

wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/007/SRR4420297/SRR4420297\_2.fastq.gz

wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/008/SRR4420298/SRR4420298\_1.fastq.gz

wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/008/SRR4420298/SRR4420298\_2.fastq.gz



cd ..



cd reference



\# Download reference genome

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF\_000001735.4\_TAIR10.1/GCF\_000001735.4\_TAIR10.1\_genomic.fna.gz



\# Download annotation file (GFF format)

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF\_000001735.4\_TAIR10.1/GCF\_000001735.4\_TAIR10.1\_genomic.gff.gz



\# Decompress files

gunzip GCF\_000001735.4\_TAIR10.1\_genomic.fna.gz

gunzip GCF\_000001735.4\_TAIR10.1\_genomic.gff.gz



cd ..



or



cd reference



\# Download GTF directly from Ensembl Plants

wget http://ftp.ensemblgenomes.org/pub/plants/release-57/gtf/arabidopsis\_thaliana/Arabidopsis\_thaliana.TAIR10.57.gtf.gz

gunzip Arabidopsis\_thaliana.TAIR10.57.gtf.gz



\# Rename for consistency (optional)

mv Arabidopsis\_thaliana.TAIR10.57.gtf GCF\_000001735.4\_TAIR10.1\_genomic.gtf



cd ..



cd raw\_data



\# Run FastQC on all samples

fastqc \*.fastq.gz -o ../qc/pre\_trim/ -t 4



\# Generate summary report with MultiQC

cd ../qc/pre\_trim

multiqc .

cd ../..



cd raw\_data



\# Trimmomatic for paired-end reads

\# Adjust adapter file path as needed (e.g., TruSeq3-PE.fa)



\# Control sample 1

trimmomatic PE -threads 4 \\

&nbsp; SRR4420293\_1.fastq.gz SRR4420293\_2.fastq.gz \\

&nbsp; ../trimmed/SRR4420293\_1\_paired.fq.gz ../trimmed/SRR4420293\_1\_unpaired.fq.gz \\

&nbsp; ../trimmed/SRR4420293\_2\_paired.fq.gz ../trimmed/SRR4420293\_2\_unpaired.fq.gz \\

&nbsp; ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36



\# Control sample 2

trimmomatic PE -threads 4 \\

&nbsp; SRR4420294\_1.fastq.gz SRR4420294\_2.fastq.gz \\

&nbsp; ../trimmed/SRR4420294\_1\_paired.fq.gz ../trimmed/SRR4420294\_1\_unpaired.fq.gz \\

&nbsp; ../trimmed/SRR4420294\_2\_paired.fq.gz ../trimmed/SRR4420294\_2\_unpaired.fq.gz \\

&nbsp; ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36



\# Test sample 1

trimmomatic PE -threads 4 \\

&nbsp; SRR4420297\_1.fastq.gz SRR4420297\_2.fastq.gz \\

&nbsp; ../trimmed/SRR4420297\_1\_paired.fq.gz ../trimmed/SRR4420297\_1\_unpaired.fq.gz \\

&nbsp; ../trimmed/SRR4420297\_2\_paired.fq.gz ../trimmed/SRR4420297\_2\_unpaired.fq.gz \\

&nbsp; ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36



\# Test sample 2

trimmomatic PE -threads 4 \\

&nbsp; SRR4420298\_1.fastq.gz SRR4420298\_2.fastq.gz \\

&nbsp; ../trimmed/SRR4420298\_1\_paired.fq.gz ../trimmed/SRR4420298\_1\_unpaired.fq.gz \\

&nbsp; ../trimmed/SRR4420298\_2\_paired.fq.gz ../trimmed/SRR4420298\_2\_unpaired.fq.gz \\

&nbsp; ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36



cd ..





\# Run FastQC on trimmed reads

fastqc trimmed/\*\_paired.fq.gz -o qc/post\_trim/ -t 4



\# Generate MultiQC report

cd qc/post\_trim

multiqc .

cd ../..





cd reference



\# Step 7.1: Extract splice sites from GTF file

echo "Extracting splice sites..."

hisat2\_extract\_splice\_sites.py GCF\_000001735.4\_TAIR10.1\_genomic.gtf > splice\_sites.txt



\# Verify splice sites were extracted

echo "Number of splice sites extracted:"

wc -l splice\_sites.txt



\# Step 7.2: Extract exons from GTF file

echo "Extracting exons..."

hisat2\_extract\_exons.py GCF\_000001735.4\_TAIR10.1\_genomic.gtf > exons.txt



\# Verify exons were extracted

echo "Number of exons extracted:"

wc -l exons.txt



\# Step 7.3: Build HISAT2 index WITH splice site and exon information

echo "Building HISAT2 index with splice-awareness..."

hisat2-build -p 4 \\

&nbsp; --ss splice\_sites.txt \\

&nbsp; --exon exons.txt \\

&nbsp; GCF\_000001735.4\_TAIR10.1\_genomic.fna \\

&nbsp; arabidopsis\_index



\# Verify index files were created

echo "HISAT2 index files created:"

ls -lh arabidopsis\_index\*.ht2



cd ..



cd trimmed

\# Control sample 1

hisat2 -p 4 -x ../reference/arabidopsis\_index \\

&nbsp; --dta \\

&nbsp; --summary-file ../alignment/SRR4420293\_summary.txt \\

&nbsp; -1 SRR4420293\_1\_paired.fq.gz -2 SRR4420293\_2\_paired.fq.gz \\

&nbsp; -S ../alignment/SRR4420293.sam

\# Control sample 2

hisat2 -p 4 -x ../reference/arabidopsis\_index \\

&nbsp; --dta \\

&nbsp; --summary-file ../alignment/SRR4420294\_summary.txt \\

&nbsp; -1 SRR4420294\_1\_paired.fq.gz -2 SRR4420294\_2\_paired.fq.gz \\

&nbsp; -S ../alignment/SRR4420294.sam

\# Test sample 1

hisat2 -p 4 -x ../reference/arabidopsis\_index \\

&nbsp; --dta \\

&nbsp; --summary-file ../alignment/SRR4420297\_summary.txt \\

&nbsp; -1 SRR4420297\_1\_paired.fq.gz -2 SRR4420297\_2\_paired.fq.gz \\

&nbsp; -S ../alignment/SRR4420297.sam

\# Test sample 2

hisat2 -p 4 -x ../reference/arabidopsis\_index \\

&nbsp; --dta \\

&nbsp; --summary-file ../alignment/SRR4420298\_summary.txt \\

&nbsp; -1 SRR4420298\_1\_paired.fq.gz -2 SRR4420298\_2\_paired.fq.gz \\

&nbsp; -S ../alignment/SRR4420298.sam

cd ..





cd alignment



\# Process each sample

for sample in SRR4420293 SRR4420294 SRR4420297 SRR4420298

do

&nbsp; echo "Processing $sample..."

&nbsp; 

&nbsp; # Convert SAM to BAM

&nbsp; samtools view -@ 4 -bS ${sample}.sam > ${sample}.bam

&nbsp; 

&nbsp; # Sort BAM file

&nbsp; samtools sort -@ 4 ${sample}.bam -o ${sample}\_sorted.bam

&nbsp; 

&nbsp; # Index BAM file

&nbsp; samtools index ${sample}\_sorted.bam

&nbsp; 

&nbsp; # Remove unsorted BAM and SAM to save space

&nbsp; rm ${sample}.bam ${sample}.sam

&nbsp; 

&nbsp; echo "$sample processing complete"

done



cd ..



cd alignment



\# Generate alignment statistics

for sample in SRR4420293 SRR4420294 SRR4420297 SRR4420298

do

&nbsp; echo "Generating statistics for $sample..."

&nbsp; samtools flagstat ${sample}\_sorted.bam > ${sample}\_stats.txt

&nbsp; samtools idxstats ${sample}\_sorted.bam > ${sample}\_idxstats.txt

done



\# Run MultiQC on alignment results

multiqc .



cd ..



cd alignment



\# Count reads using featureCounts with the GTF file

featureCounts -T 4 -p -t exon -g gene\_id \\

&nbsp; -a ../reference/GCF\_000001735.4\_TAIR10.1\_genomic.gtf \\

&nbsp; -o ../counts/gene\_counts.txt \\

&nbsp; SRR4420293\_sorted.bam \\

&nbsp; SRR4420294\_sorted.bam \\

&nbsp; SRR4420297\_sorted.bam \\

&nbsp; SRR4420298\_sorted.bam



\# Generate summary statistics

cd ../counts

echo "Gene counting complete. Summary:"

head -n 50 gene\_counts.txt.summary



cd ..



**# Load required libraries**

**library(DESeq2)**

**library(ggplot2)**

**library(pheatmap)**

**library(EnhancedVolcano)**



**# Set working directory**

**setwd("rnaseq\_analysis/counts")**



**# Read count data**

**countData <- read.table("gene\_counts.txt", header=TRUE, row.names=1, skip=1)**



**# Remove first 5 columns (Chr, Start, End, Strand, Length)**

**countData <- countData\[, 6:ncol(countData)]**



**# Simplify column names**

**colnames(countData) <- c("Control1", "Control2", "Test1", "Test2")**



**# Create sample information**

**colData <- data.frame(**

  **condition = factor(c("Control", "Control", "Test", "Test")),**

  **row.names = colnames(countData)**

**)**



**# Create DESeq2 dataset**

**dds <- DESeqDataSetFromMatrix(**

  **countData = countData,**

  **colData = colData,**

  **design = ~ condition**

**)**



**# Pre-filtering: remove genes with very low counts**

**keep <- rowSums(counts(dds)) >= 10**

**dds <- dds\[keep,]**



**# Run DESeq2 analysis**

**dds <- DESeq(dds)**



**# Get results**

**res <- results(dds, contrast=c("condition", "Test", "Control"))**



**# Order by adjusted p-value**

**resOrdered <- res\[order(res$padj),]**



**# Summary of results**

**summary(res)**



**# Save results**

**write.csv(as.data.frame(resOrdered),** 

          **file="../results/differential\_expression\_results.csv")**



**# Extract significant genes (padj < 0.05, |log2FC| > 1)**

**sig\_genes <- subset(resOrdered, padj < 0.05 \& abs(log2FoldChange) > 1)**

**write.csv(as.data.frame(sig\_genes),** 

          **file="../results/significant\_genes.csv")**



**# Normalized counts for visualization**

**normalized\_counts <- counts(dds, normalized=TRUE)**

**write.csv(normalized\_counts,** 

          **file="../results/normalized\_counts.csv")**



**# Variance stabilizing transformation for visualization**

**vsd <- vst(dds, blind=FALSE)**



**# PCA plot**

**pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)**

**percentVar <- round(100 \* attr(pcaData, "percentVar"))**

**pdf("../results/PCA\_plot.pdf", width=8, height=6)**

**ggplot(pcaData, aes(PC1, PC2, color=condition, label=name)) +**

  **geom\_point(size=3) +**

  **geom\_text(vjust=1.5) +**

  **xlab(paste0("PC1: ",percentVar,"% variance")) +**

  **ylab(paste0("PC2: ",percentVar,"% variance")) +**

  **theme\_bw()**

**dev.off()**



**# Heatmap of top 50 genes**

**top\_genes <- head(order(resOrdered$padj), 50)**

**mat <- assay(vsd)\[top\_genes, ]**

**mat <- mat - rowMeans(mat)**

**pdf("../results/heatmap\_top50.pdf", width=10, height=12)**

**pheatmap(mat,** 

         **annotation\_col = colData,**

         **cluster\_rows=TRUE,** 

         **cluster\_cols=TRUE,**

         **show\_rownames=FALSE)**

**dev.off()**



**# Volcano plot**

**pdf("../results/volcano\_plot.pdf", width=10, height=8)**

**EnhancedVolcano(res,**

    **lab = rownames(res),**

    **x = 'log2FoldChange',**

    **y = 'padj',**

    **title = 'Test vs Control',**

    **pCutoff = 0.05,**

    **FCcutoff = 1.0,**

    **pointSize = 2.0,**

    **labSize = 4.0)**

**dev.off()**



**# MA plot**

**pdf("../results/MA\_plot.pdf", width=8, height=6)**

**plotMA(res, ylim=c(-5,5))**

**dev.off()**



**# Dispersion plot**

**pdf("../results/dispersion\_plot.pdf", width=8, height=6)**

**plotDispEsts(dds)**

**dev.off()**



**print("Analysis complete! Check the results folder.")**







