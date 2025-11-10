# Complete RNA-seq Data Analysis Pipeline
## With GFF to GTF Conversion and Splice Site Extraction

### Overview
- **Organism**: *Arabidopsis thaliana* (TAIR10.1)
- **Design**: 2 control samples vs 2 test samples
- **Sequencing**: Paired-end reads
- **Analysis Goal**: Differential gene expression
- **Key Fix**: Proper GFF to GTF conversion for splice site and exon extraction

---

## Step 1: Set Up Working Environment

```bash
# Create project directory structure
mkdir -p rnaseq_analysis/{raw_data,qc/{pre_trim,post_trim},trimmed,reference,alignment,counts,results}
cd rnaseq_analysis

# Install required tools (if not already installed)
# FastQC, Trimmomatic, HISAT2, SAMtools, featureCounts, DESeq2
# Additionally: gffread for GFF to GTF conversion
```

---

## Step 2: Download Raw Data

```bash
cd raw_data

# Control samples
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/003/SRR4420293/SRR4420293_1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/003/SRR4420293/SRR4420293_2.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/004/SRR4420294/SRR4420294_1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/004/SRR4420294/SRR4420294_2.fastq.gz

# Test/Treatment samples
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/007/SRR4420297/SRR4420297_1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/007/SRR4420297/SRR4420297_2.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/008/SRR4420298/SRR4420298_1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/008/SRR4420298/SRR4420298_2.fastq.gz

cd ..
```

---

## Step 3: Download Reference Genome and Annotation

```bash
cd reference

# Download reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz

# Download annotation file (GFF format)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff.gz

# Decompress files
gunzip GCF_000001735.4_TAIR10.1_genomic.fna.gz
gunzip GCF_000001735.4_TAIR10.1_genomic.gff.gz

cd ..
```

---

## Step 3.5: **NEW - Convert GFF to GTF Format**

**This is the critical step that fixes the issue your colleagues faced!**

### Method 1: Using gffread (Recommended)

```bash
cd reference

# Install gffread if not available
# conda install -c bioconda gffread
# OR: brew install gffread
# OR: compile from https://github.com/gpertea/gffread

# Convert GFF3 to GTF format
gffread GCF_000001735.4_TAIR10.1_genomic.gff -T -o GCF_000001735.4_TAIR10.1_genomic.gtf

# Verify the conversion
echo "Lines in original GFF file:"
wc -l GCF_000001735.4_TAIR10.1_genomic.gff
echo "Lines in converted GTF file:"
wc -l GCF_000001735.4_TAIR10.1_genomic.gtf

# Check for problematic characters
echo "Checking for '?' in GFF file:"
grep -c '?' GCF_000001735.4_TAIR10.1_genomic.gff || echo "No '?' found"
echo "Checking for '?' in GTF file:"
grep -c '?' GCF_000001735.4_TAIR10.1_genomic.gtf || echo "No '?' found"

cd ..
```

### Method 2: Alternative - Download Pre-converted GTF from Ensembl

```bash
cd reference

# Download GTF directly from Ensembl Plants
wget http://ftp.ensemblgenomes.org/pub/plants/release-57/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.57.gtf.gz
gunzip Arabidopsis_thaliana.TAIR10.57.gtf.gz

# Rename for consistency (optional)
mv Arabidopsis_thaliana.TAIR10.57.gtf GCF_000001735.4_TAIR10.1_genomic.gtf

cd ..
```

### Method 3: Python Script for Manual Conversion

If gffread is not available, create this Python script:

```bash
cd reference
cat > convert_gff_to_gtf.py << 'EOF'
#!/usr/bin/env python3
"""
Convert GFF3 to GTF format and remove problematic lines
This handles the issue with '?' characters and malformed lines
"""

def clean_and_convert_gff_to_gtf(input_gff, output_gtf):
    with open(input_gff, 'r') as infile, open(output_gtf, 'w') as outfile:
        line_count = 0
        skipped_count = 0
        
        for line in infile:
            # Skip comment lines
            if line.startswith('#'):
                continue
            
            # Skip FASTA section
            if line.startswith('>'):
                break
            
            # Skip empty lines
            if not line.strip():
                continue
            
            # Skip lines with '?' which cause parsing errors
            if '?' in line:
                skipped_count += 1
                continue
            
            fields = line.strip().split('\t')
            
            # GFF must have exactly 9 fields
            if len(fields) != 9:
                skipped_count += 1
                continue
            
            # Extract fields
            seqname, source, feature, start, end, score, strand, frame, attributes = fields
            
            # Filter for relevant features for RNA-seq analysis
            # Keep gene, transcript, exon, CDS for proper annotation
            if feature in ['gene', 'transcript', 'mRNA', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR']:
                # Write to GTF file (GTF is same format as GFF for these purposes)
                outfile.write('\t'.join(fields) + '\n')
                line_count += 1
        
        print(f"Conversion complete!")
        print(f"Lines written: {line_count}")
        print(f"Lines skipped: {skipped_count}")

if __name__ == "__main__":
    input_file = "GCF_000001735.4_TAIR10.1_genomic.gff"
    output_file = "GCF_000001735.4_TAIR10.1_genomic.gtf"
    clean_and_convert_gff_to_gtf(input_file, output_file)
EOF

# Make executable and run
chmod +x convert_gff_to_gtf.py
python3 convert_gff_to_gtf.py

cd ..
```

---

## Step 4: Quality Control (Pre-trimming)

```bash
cd raw_data

# Run FastQC on all samples
fastqc *.fastq.gz -o ../qc/pre_trim/ -t 4

# Generate summary report with MultiQC
cd ../qc/pre_trim
multiqc .
cd ../..
```

**What to check:**
- Per base sequence quality
- Adapter content
- Overrepresented sequences
- GC content distribution

---

## Step 5: Read Trimming and Quality Filtering

```bash
cd raw_data

# Trimmomatic for paired-end reads
# Adjust adapter file path as needed (e.g., TruSeq3-PE.fa)

# Control sample 1
trimmomatic PE -threads 4 \
  SRR4420293_1.fastq.gz SRR4420293_2.fastq.gz \
  ../trimmed/SRR4420293_1_paired.fq.gz ../trimmed/SRR4420293_1_unpaired.fq.gz \
  ../trimmed/SRR4420293_2_paired.fq.gz ../trimmed/SRR4420293_2_unpaired.fq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Control sample 2
trimmomatic PE -threads 4 \
  SRR4420294_1.fastq.gz SRR4420294_2.fastq.gz \
  ../trimmed/SRR4420294_1_paired.fq.gz ../trimmed/SRR4420294_1_unpaired.fq.gz \
  ../trimmed/SRR4420294_2_paired.fq.gz ../trimmed/SRR4420294_2_unpaired.fq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Test sample 1
trimmomatic PE -threads 4 \
  SRR4420297_1.fastq.gz SRR4420297_2.fastq.gz \
  ../trimmed/SRR4420297_1_paired.fq.gz ../trimmed/SRR4420297_1_unpaired.fq.gz \
  ../trimmed/SRR4420297_2_paired.fq.gz ../trimmed/SRR4420297_2_unpaired.fq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Test sample 2
trimmomatic PE -threads 4 \
  SRR4420298_1.fastq.gz SRR4420298_2.fastq.gz \
  ../trimmed/SRR4420298_1_paired.fq.gz ../trimmed/SRR4420298_1_unpaired.fq.gz \
  ../trimmed/SRR4420298_2_paired.fq.gz ../trimmed/SRR4420298_2_unpaired.fq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

cd ..
```

---

## Step 6: Post-trimming Quality Control

```bash
# Run FastQC on trimmed reads
fastqc trimmed/*_paired.fq.gz -o qc/post_trim/ -t 4

# Generate MultiQC report
cd qc/post_trim
multiqc .
cd ../..
```

---

## Step 7: **REVISED - Extract Splice Sites and Exons, then Build HISAT2 Index**

**This step now uses the converted GTF file to properly extract splice sites and exons!**

```bash
cd reference

# Step 7.1: Extract splice sites from GTF file
echo "Extracting splice sites..."
hisat2_extract_splice_sites.py GCF_000001735.4_TAIR10.1_genomic.gtf > splice_sites.txt

# Verify splice sites were extracted
echo "Number of splice sites extracted:"
wc -l splice_sites.txt

# Step 7.2: Extract exons from GTF file
echo "Extracting exons..."
hisat2_extract_exons.py GCF_000001735.4_TAIR10.1_genomic.gtf > exons.txt

# Verify exons were extracted
echo "Number of exons extracted:"
wc -l exons.txt

# Step 7.3: Build HISAT2 index WITH splice site and exon information
echo "Building HISAT2 index with splice-awareness..."
hisat2-build -p 4 \
  --ss splice_sites.txt \
  --exon exons.txt \
  GCF_000001735.4_TAIR10.1_genomic.fna \
  arabidopsis_index

# Verify index files were created
echo "HISAT2 index files created:"
ls -lh arabidopsis_index*.ht2

cd ..
```

**What this does:**
- `hisat2_extract_splice_sites.py`: Extracts known splice junction sites from the GTF file
- `hisat2_extract_exons.py`: Extracts exon boundaries from the GTF file
- `--ss` and `--exon` flags: Incorporate this information into the HISAT2 index
- **Result**: Better alignment accuracy, especially for reads spanning splice junctions

---

## Step 8: Alignment to Reference Genome

```bash
cd trimmed

# Control sample 1
hisat2 -p 4 -x ../reference/arabidopsis_index \
  --dta \
  --summary-file ../alignment/SRR4420293_summary.txt \
  -1 SRR4420293_1_paired.fq.gz -2 SRR4420293_2_paired.fq.gz \
  -S ../alignment/SRR4420293.sam

# Control sample 2
hisat2 -p 4 -x ../reference/arabidopsis_index \
  --dta \
  --summary-file ../alignment/SRR4420294_summary.txt \
  -1 SRR4420294_1_paired.fq.gz -2 SRR4420294_2_paired.fq.gz \
  -S ../alignment/SRR4420294.sam

# Test sample 1
hisat2 -p 4 -x ../reference/arabidopsis_index \
  --dta \
  --summary-file ../alignment/SRR4420297_summary.txt \
  -1 SRR4420297_1_paired.fq.gz -2 SRR4420297_2_paired.fq.gz \
  -S ../alignment/SRR4420297.sam

# Test sample 2
hisat2 -p 4 -x ../reference/arabidopsis_index \
  --dta \
  --summary-file ../alignment/SRR4420298_summary.txt \
  -1 SRR4420298_1_paired.fq.gz -2 SRR4420298_2_paired.fq.gz \
  -S ../alignment/SRR4420298.sam

cd ..
```

**Parameters explained:**
- `--dta`: Report alignments tailored for transcript assemblers (downstream transcriptome assembly)
- `--summary-file`: Write alignment summary statistics to file

---

## Step 9: SAM to BAM Conversion and Sorting

```bash
cd alignment

# Process each sample
for sample in SRR4420293 SRR4420294 SRR4420297 SRR4420298
do
  echo "Processing $sample..."
  
  # Convert SAM to BAM
  samtools view -@ 4 -bS ${sample}.sam > ${sample}.bam
  
  # Sort BAM file
  samtools sort -@ 4 ${sample}.bam -o ${sample}_sorted.bam
  
  # Index BAM file
  samtools index ${sample}_sorted.bam
  
  # Remove unsorted BAM and SAM to save space
  rm ${sample}.bam ${sample}.sam
  
  echo "$sample processing complete"
done

cd ..
```

---

## Step 10: Alignment Quality Assessment

```bash
cd alignment

# Generate alignment statistics
for sample in SRR4420293 SRR4420294 SRR4420297 SRR4420298
do
  echo "Generating statistics for $sample..."
  samtools flagstat ${sample}_sorted.bam > ${sample}_stats.txt
  samtools idxstats ${sample}_sorted.bam > ${sample}_idxstats.txt
done

# Run MultiQC on alignment results
multiqc .

cd ..
```

---

## Step 11: **REVISED - Read Quantification Using GTF File**

```bash
cd alignment

# Count reads using featureCounts with the GTF file
featureCounts -T 4 -p -t exon -g gene_id \
  -a ../reference/GCF_000001735.4_TAIR10.1_genomic.gtf \
  -o ../counts/gene_counts.txt \
  SRR4420293_sorted.bam \
  SRR4420294_sorted.bam \
  SRR4420297_sorted.bam \
  SRR4420298_sorted.bam

# Generate summary statistics
cd ../counts
echo "Gene counting complete. Summary:"
head -n 50 gene_counts.txt.summary

cd ..
```

**Parameters explained:**
- `-T 4`: Use 4 threads
- `-p`: Count fragments (paired-end mode)
- `-t exon`: Count at exon level
- `-g gene_id`: Group by gene_id attribute
- **Now using GTF file instead of GFF!**

---

## Step 12: Differential Expression Analysis (DESeq2 in R)

Create an R script `deseq2_analysis.R`:

```R
# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

# Set working directory
setwd("rnaseq_analysis/counts")

# Read count data
countData <- read.table("gene_counts.txt", header=TRUE, row.names=1, skip=1)

# Remove first 5 columns (Chr, Start, End, Strand, Length)
countData <- countData[, 6:ncol(countData)]

# Simplify column names
colnames(countData) <- c("Control1", "Control2", "Test1", "Test2")

# Create sample information
colData <- data.frame(
  condition = factor(c("Control", "Control", "Test", "Test")),
  row.names = colnames(countData)
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colData,
  design = ~ condition
)

# Pre-filtering: remove genes with very low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast=c("condition", "Test", "Control"))

# Order by adjusted p-value
resOrdered <- res[order(res$padj),]

# Summary of results
summary(res)

# Save results
write.csv(as.data.frame(resOrdered), 
          file="../results/differential_expression_results.csv")

# Extract significant genes (padj < 0.05, |log2FC| > 1)
sig_genes <- subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(sig_genes), 
          file="../results/significant_genes.csv")

# Normalized counts for visualization
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, 
          file="../results/normalized_counts.csv")

# Variance stabilizing transformation for visualization
vsd <- vst(dds, blind=FALSE)

# PCA plot
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf("../results/PCA_plot.pdf", width=8, height=6)
ggplot(pcaData, aes(PC1, PC2, color=condition, label=name)) +
  geom_point(size=3) +
  geom_text(vjust=1.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw()
dev.off()

# Heatmap of top 50 genes
top_genes <- head(order(resOrdered$padj), 50)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)
pdf("../results/heatmap_top50.pdf", width=10, height=12)
pheatmap(mat, 
         annotation_col = colData,
         cluster_rows=TRUE, 
         cluster_cols=TRUE,
         show_rownames=FALSE)
dev.off()

# Volcano plot
pdf("../results/volcano_plot.pdf", width=10, height=8)
EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    title = 'Test vs Control',
    pCutoff = 0.05,
    FCcutoff = 1.0,
    pointSize = 2.0,
    labSize = 4.0)
dev.off()

# MA plot
pdf("../results/MA_plot.pdf", width=8, height=6)
plotMA(res, ylim=c(-5,5))
dev.off()

# Dispersion plot
pdf("../results/dispersion_plot.pdf", width=8, height=6)
plotDispEsts(dds)
dev.off()

print("Analysis complete! Check the results folder.")
```

Run the R script:

```bash
Rscript deseq2_analysis.R
```

---

## Complete Automated Bash Script

Save this as `rnaseq_pipeline_fixed.sh`:

```bash
#!/bin/bash
# RNA-seq Analysis Pipeline with GFF to GTF Conversion
# Fixed version addressing splice site extraction issues
# Date: $(date)

set -e  # Exit on error

echo "========================================="
echo "RNA-seq Analysis Pipeline - Fixed Version"
echo "========================================="

# Step 1: Setup
echo "[Step 1] Setting up directories..."
mkdir -p rnaseq_analysis/{raw_data,qc/{pre_trim,post_trim},trimmed,reference,alignment,counts,results}
cd rnaseq_analysis

# Step 2: Download data
echo "[Step 2] Downloading raw data..."
cd raw_data
wget -q http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/003/SRR4420293/SRR4420293_1.fastq.gz &
wget -q http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/003/SRR4420293/SRR4420293_2.fastq.gz &
wget -q http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/004/SRR4420294/SRR4420294_1.fastq.gz &
wget -q http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/004/SRR4420294/SRR4420294_2.fastq.gz &
wget -q http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/007/SRR4420297/SRR4420297_1.fastq.gz &
wget -q http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/007/SRR4420297/SRR4420297_2.fastq.gz &
wget -q http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/008/SRR4420298/SRR4420298_1.fastq.gz &
wget -q http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/008/SRR4420298/SRR4420298_2.fastq.gz &
wait
cd ..

# Step 3: Download reference
echo "[Step 3] Downloading reference genome and annotation..."
cd reference
wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff.gz
gunzip *.gz

# Step 3.5: CRITICAL FIX - Convert GFF to GTF
echo "[Step 3.5] Converting GFF to GTF format..."
if command -v gffread &> /dev/null; then
    echo "Using gffread for conversion..."
    gffread GCF_000001735.4_TAIR10.1_genomic.gff -T -o GCF_000001735.4_TAIR10.1_genomic.gtf
else
    echo "gffread not found. Downloading pre-converted GTF from Ensembl..."
    wget -q http://ftp.ensemblgenomes.org/pub/plants/release-57/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.57.gtf.gz
    gunzip Arabidopsis_thaliana.TAIR10.57.gtf.gz
    mv Arabidopsis_thaliana.TAIR10.57.gtf GCF_000001735.4_TAIR10.1_genomic.gtf
fi
echo "GTF file created successfully!"
cd ..

# Step 4: QC pre-trim
echo "[Step 4] Running pre-trim QC..."
fastqc raw_data/*.fastq.gz -o qc/pre_trim/ -t 4 -q
cd qc/pre_trim && multiqc . && cd ../..

# Step 5: Trimming
echo "[Step 5] Trimming reads..."
for sample in SRR4420293 SRR4420294 SRR4420297 SRR4420298
do
  trimmomatic PE -threads 4 \
    raw_data/${sample}_1.fastq.gz raw_data/${sample}_2.fastq.gz \
    trimmed/${sample}_1_paired.fq.gz trimmed/${sample}_1_unpaired.fq.gz \
    trimmed/${sample}_2_paired.fq.gz trimmed/${sample}_2_unpaired.fq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

# Step 6: QC post-trim
echo "[Step 6] Running post-trim QC..."
fastqc trimmed/*_paired.fq.gz -o qc/post_trim/ -t 4 -q
cd qc/post_trim && multiqc . && cd ../..

# Step 7: Extract splice sites and exons, then build index
echo "[Step 7] Extracting splice sites and exons from GTF..."
cd reference
hisat2_extract_splice_sites.py GCF_000001735.4_TAIR10.1_genomic.gtf > splice_sites.txt
hisat2_extract_exons.py GCF_000001735.4_TAIR10.1_genomic.gtf > exons.txt

echo "Building HISAT2 index with splice-awareness..."
hisat2-build -p 4 \
  --ss splice_sites.txt \
  --exon exons.txt \
  GCF_000001735.4_TAIR10.1_genomic.fna \
  arabidopsis_index
cd ..

# Step 8: Alignment
echo "[Step 8] Aligning reads..."
for sample in SRR4420293 SRR4420294 SRR4420297 SRR4420298
do
  hisat2 -p 4 -x reference/arabidopsis_index \
    --dta \
    --summary-file alignment/${sample}_summary.txt \
    -1 trimmed/${sample}_1_paired.fq.gz \
    -2 trimmed/${sample}_2_paired.fq.gz \
    -S alignment/${sample}.sam
done

# Step 9: SAM to BAM
echo "[Step 9] Converting to BAM and sorting..."
cd alignment
for sample in SRR4420293 SRR4420294 SRR4420297 SRR4420298
do
  samtools view -@ 4 -bS ${sample}.sam | samtools sort -@ 4 -o ${sample}_sorted.bam
  samtools index ${sample}_sorted.bam
  rm ${sample}.sam
done
cd ..

# Step 10: Alignment stats
echo "[Step 10] Generating alignment statistics..."
cd alignment
for sample in SRR4420293 SRR4420294 SRR4420297 SRR4420298
do
  samtools flagstat ${sample}_sorted.bam > ${sample}_stats.txt
done
multiqc .
cd ..

# Step 11: Count reads using GTF file
echo "[Step 11] Counting reads..."
featureCounts -T 4 -p -t exon -g gene_id \
  -a reference/GCF_000001735.4_TAIR10.1_genomic.gtf \
  -o counts/gene_counts.txt \
  alignment/*_sorted.bam

echo "========================================="
echo "Pipeline complete! Run DESeq2 analysis in R."
echo "========================================="
```

Make it executable:
```bash
chmod +x rnaseq_pipeline_fixed.sh
./rnaseq_pipeline_fixed.sh
```

---

## Summary of Key Changes

### The Problem
- **Original issue**: GFF file contained '?' characters and formatting that caused `hisat2_extract_splice_sites.py` and `hisat2_extract_exons.py` to fail
- The hisat Python scripts require strict GTF format

### The Solution
1. **Convert GFF to GTF** using `gffread` or download pre-converted GTF from Ensembl
2. **Extract splice sites and exons** from the clean GTF file
3. **Build HISAT2 index** with splice-awareness using `--ss` and `--exon` flags
4. **Use GTF file** for featureCounts quantification

### Benefits of This Approach
- ✅ Eliminates parsing errors from '?' characters
- ✅ Improves alignment accuracy at splice junctions
- ✅ Better handling of intron-spanning reads
- ✅ More accurate gene expression quantification
- ✅ Compatible with all downstream analysis tools

---

## Troubleshooting

### If splice site extraction still fails:
```bash
# Check GTF file format
head -n 20 reference/GCF_000001735.4_TAIR10.1_genomic.gtf

# Check for problematic characters
grep '?' reference/GCF_000001735.4_TAIR10.1_genomic.gtf

# Manually test extraction
cd reference
python3 $(which hisat2_extract_splice_sites.py) GCF_000001735.4_TAIR10.1_genomic.gtf | head
```

### If gffread is not available:
Install via conda:
```bash
conda install -c bioconda gffread
```

Or download pre-converted GTF directly as shown in Method 2.

---

## Expected Output Files

After completing the pipeline, you should have:
- ✅ `reference/GCF_000001735.4_TAIR10.1_genomic.gtf` - Converted GTF file
- ✅ `reference/splice_sites.txt` - Extracted splice sites
- ✅ `reference/exons.txt` - Extracted exons
- ✅ `reference/arabidopsis_index*.ht2` - HISAT2 index files (8 files)
- ✅ `alignment/*_sorted.bam` - Aligned and sorted BAM files
- ✅ `counts/gene_counts.txt` - Gene expression counts
- ✅ `results/differential_expression_results.csv` - DESeq2 results

---

## Contact and Support

For issues with:
- GFF/GTF conversion: Check gffread documentation
- HISAT2 indexing: Verify GTF format compliance
- Downstream analysis: Consult DESeq2 vignettes

**This pipeline resolves the issue faced by your colleagues and should work correctly!**
