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


chmod +x rnaseq_pipeline_fixed.sh
./rnaseq_pipeline_fixed.sh
