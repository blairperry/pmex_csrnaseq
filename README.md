
# Analysis of transcription in H2S-adapted fish (*Poecilia mexicana*)

This repo contains code for all data processing and analysis used in our study of transcription and gene expression in H2S-adapted populations of *P. mexicana*.

For questions, please contact blair.perry(at)wsu.edu.

## Data Availability
- Newly-generated csRNA-seq data
	- Raw data available at NCBI BioProject: PRJNA743555
- Previously published mRNA-seq data ([Kelley et al. 2016](http://mbe.oxfordjournals.org/content/early/2016/03/01/molbev.msw020))
	- Raw data available at NCBI BioProject: [PRJNA290391](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA290391/)

## Contents
1. [cRNA-seq trimming, mapping, peak calling, and peak annotation](#1-crna-seq-qc-mapping-peak-calling-and-annotation)
2. [mRNA-seq trimming and mapping ](#2-mrna-seq-trimming-and-mapping)
3. [Analyses of differential transcript initiation (csRNA)](#3-analyses-of-differential-transcript-initiation-csrna)
4. [Analyses of differential gene expression (mRNA)](#4-analyses-of-differential-gene-expression-mrna)
5. [Joint analyses of csRNA and mRNA data](#5-joint-analyses-of-csrna-and-mrna-data)
	a. [Characterization of transcript initiation annotations](#a-characterization-of-transcript-initiation-annotations)
	b. [Identification of putative enhancer regions](#b-identification-and-analysis-of-putative-enhancer-regions)
6. TF binding site enrichment analyses


---
## 1. cRNA-seq QC, mapping, peak calling and annotation

These analyses were run on WSU's HPC ([Kamiak](https://hpc.wsu.edu/)). For generalizability, simplified commands are presented here rather than the specific SLURM scripts used to run these commands on Kamiak. 

Note: The variable `$sample` should be replaced with the full sample name.
### Trim with homerTools *trim*
```bash
homerTools trim \
  -3 AGATCGGAAGAGCACACGTCT \
  -mis 2 \
  -minMatchLength 4 \
  -min 20 \
  -q 20 \
  $sample.fastq.gz
```

### Map csRNA-seq with STAR
Reference genome and annotations in GFF format were downloaded from NCBI (Accession: [GCF_001443325.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001443325.1/)).

```bash
# Index genome for use with STAR
STAR \
  --runThreadN 7 \
  --runMode genomeGenerate \
  --genomeDir genomeDir \
  --genomeFastaFiles GCF_001443325.1_P_mexicana-1.0_genomic.fna \
  --sjdbGTFfile GCF_001443325.1_P_mexicana-1.0_genomic.gff \
  --sjdbOverhang 83 \
  --genomeSAindexNbases 13

# Map Reads
STAR \
	--genomeDir genomeDir \
	--runThreadN 24 \
	--readFilesIn $sample.fastq.gz.trimmed \
	--outFileNamePrefix $sample. \
	--outSAMstrandField intronMotif \
	--outMultimapperOrder Random \
	--outSAMmultNmax 1 \
	--outFilterMultimapNmax 10000 \
	--limitOutSAMoneReadBytes 10000000
```

### Make tag directories
Create tag directories for csRNA-seq, input, and mRNA-seq data (see below for details on mRNA-seq processing)

```bash
makeTagDirectory \
  $sample\_tagDir \
  $sample.Aligned.out.sam \
  -genome $reference_genome \
  -checkGC \
  -unique \
  -fragLength 30 \
  -single
```

### Identify transcription initiation peaks
Peaks are clusters of transcription initiation, called using the "full approach" of the HOMER pipeline. This incorporates mRNA tag directories into the csRNA peak calling to help eliminate false positives like exomic contaminants.

```bash
# Convert GFF to GTF format
gffread \
	-T GCF_001443325.1_P_mexicana-1.0_genomic.gff \
	-o GCF_001443325.1_P_mexicana-1.0_genomic.gtf

findcsRNATSS.pl \
  $sample\_csRNA_tagDir  \
  -o $sample \
  -i $sample\_input_tagDir \
  -rna $sample\_rna_tagDir \
  -gtf $annotations \
  -genome GCF_001443325.1_P_mexicana-1.0_genomic.fna \
  -ntagThreshold 7
```

### Merge peaks
Produce a non-redundant list of peaks from all samples; 12 samples in this study.

```bash
mergePeaks \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	-strand \
	> merged.tss.txt
```

### Annotate peaks
Quantify csRNAseq read counts and annotations for merged peaks. 

```bash
annotatePeaks.pl \
	merged.tss.txt \
	GCF_001443325.1_P_mexicana-1.0_genomic.fna \
	-gtf GCF_001443325.1_P_mexicana-1.0_genomic.gtf \
	-strand + \
	-fragLength 1 \
	-raw \
	-d $sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	> rawcounts.txt
```

---

## 2. mRNA-seq trimming and mapping 
The variable $sample should be replaced with the full sample name.
### Trim with TrimGalore
```bash
# Trim adapters
trim_galore \
  --output_dir /path/to/outdir \
  --quality 0 \
  --fastqc \
  --fastqc_args "--nogroup --noextract --outdir /path/to/outdir" \
  -path_to_cutadapt /path/to/cutadapt \
  --illumina \
  --stringency 6 \
  --clip_R1 11 \
  --clip_R2 11 \
  --paired \
  --gzip ${sample}_1.fq.gz ${sample}_2.fq.gz

# Quality trim
trim_galore \
  --output_dir /path/to/outdir \
  --quality 24 \
  --fastqc \
  --fastqc_args "--nogroup --noextract --outdir /path/to/outdir" \
  -path_to_cutadapt /path/to/cutadapt \
  --illumina \
  --stringency 6 \
  --length 50 \
  --paired \
  --gzip ${sample}_1_val_1.fq.gz ${sample}_2_val_2.fq.gz
```


### Map RNA-seq with HISAT2

```bash

# Index the genome
hisat2-build \
  -f GCF_001443325.1_P_mexicana-1.0_genomic_with_mito_sequence.fa \
  GCF_001443325.1_P_mexicana-1.0_genomic_with_mito_sequence

# Map 
hisat2 \
  -q \
  --threads 4 \
  -k 5 \
  --downstream-transcriptome-assembly \
  --fr \
  --rg SM:${sample} \
  --rg ID:${rgid} \
  --rg PL:ILLUMINA \
  -x GCF_001443325.1_P_mexicana-1.0_genomic_with_mito_sequence \
  -1 ${sample}_1_val_1_val_1.fq.gz \
  -2 ${sample}_2_val_2_val_2.fq.gz \
  -S ${sample}_HISat_out.sam \
  --summary-file ${sample}_statistics.txt
```

### Convert, sort, and merge SAM/BAM files

```bash

# Convert to BAM
samtools view \
  -bSh ${sample}_HISat_out.sam > \
  ${sample}_HISat_out.bam
  
# Sort by coordinate
samtools sort \
  -@ 8 \
  ${sample}_HISat_out.bam \
  -o ${sample}_sorted.bam
  
# Merge BAM files from same individual
samtools merge \
  -r \
  ${sample}_merged_sorted.bam \
  ${sample}_sorted.bam \
  ${sample}_sorted.bam

```

### Build gene count matrices

Note: modified prepDE.py script ([prepDE-2019-10-01_edited.py](https://github.com/blairperry/pmex_csrnaseq/blob/main/prepDE-2019-10-01_edited.py)) is provided in this repository. 

```bash
# Generate GTF files for each sample using StringTie and modified version of the associated prepDE.py script

mkdir -p ballgown

stringtie \
  ${sample}_sorted.bam \
  -o ballgown/${sample}/${sample}.gtf \
  -p 8 \
  -G GCF_001443325.1_P_mexicana-1.0_genomic_with_mito_sequence-edited.gff \
  -B \
  -e

# Generate gene count matrix
grep -v "STRG" gene_count_matrix.csv > gene_count_matrix_no_STRG.csv

```

---
## 3. Analyses of differential transcript initiation (csRNA) 

The following R script contains code used to:
- Filter csRNA-seq count data 
- Perform analysis of differential transcript initiation between sulfidic and non-sulfidic populations using edgeR
- Plot heatmaps of csRNA-seq peaks for Figure 2. 
- Assess differential initiation in candidate H2S detox genes 

Link to Rscript: [3_csRNA_DiffInitiation_01.15.24.R](https://github.com/blairperry/pmex_csrnaseq/blob/main/analysis/3_csRNA_DiffInitiation_01.15.24.R)

---
## 4. Analyses of differential gene expression (mRNA)

The following R script contains code used to:
- Filter mRNA-seq count data 
- Perform analysis of differential gene expression between sulfidic and non-sulfidic populations using edgeR

Link to Rscript: [4_mRNA_DiffExpression_01.15.24.R](https://github.com/blairperry/pmex_csrnaseq/blob/main/analysis/4_mRNA_DiffExpression_01.15.24.R)

---
## 5. Joint analyses of csRNA and mRNA data
### a. Characterization of transcript initiation annotations
The following R script contains code used to:
- Summarize annotations of differentially initiated peaks (i.e., whether peaks are located in promoters, exons, etc.)
- Identify differentially initiated peaks (csRNA) associated with differentially expressed genes (mRNA)
- Test whether specific subsets of differentially initiated peaks (e.g., DI peaks associated with DE genes, DI peaks associated with non-DE genes) exhibit unique distributions of peak annotations
- Generate barplots and alluvial plots for Figure 3. 

Link to Rscript: [5a_peakAnnotationExploration_01.15.24.R](https://github.com/blairperry/pmex_csrnaseq/blob/main/analysis/5a_peakAnnotationExploration_01.15.24.R)

### b. Identification and analysis of putative enhancer regions

The following R script contains code used to:
- Parse and summarize csRNA peak annotations
- Identify putative enhancer regions defined as csRNAseq peaks that are distal from protein coding genes and comprised of unstable transcripts

Link to Rscript: [5b_1_putEnhIdentification_01.16.24.R](https://github.com/blairperry/pmex_csrnaseq/blob/main/analysis/5b_1_putEnhIdentification_01.16.24.R)

Using the bedfile of all putative enhancer regions identified above, the following commands were used to expand each region by 500 kbp in each direction and intersect with a bedfile of gene TSS positions. Genes found to be located within 500 kbp of a putative enhancer are considered potential target genes that are evaluated further below. 

```bash
# Sort bedfile of putative enhancer regions
bedtools sort -i putativeEnhancerPeaks_pmex_01.12.24.bed > putativeEnhancerPeaks_pmex_01.12.24.sort.bed

# Expand each putative enhancer region by 500 kbp in both directions
bedtools slop -b 500000 -i putativeEnhancerPeaks_pmex_01.12.24.sort.bed -g ../../data/reference/chrNameLength.txt > putativeEnhancerPeaks_pmex_01.12.24.500kbSlop.bed

# Intersect with a bed file of gene TSS positions 
bedtools intersect -a putativeEnhancerPeaks_pmex_01.12.24.500kbSlop.bed -b ../../data/reference/tss.bed -wa -wb > putativeEnhancerPeaks_pmex_01.12.24.500kbSlop.tssIntersect.txt
```

The following R script contains code used to:
- Identify putative enhancer regions with evidence of differential initiation (i.e., differential activity) in csRNA-seq edgeR analyses described above.
- Identify putative enhancers and their putative target genes with positively correlated log2 fold-changes (i.e., shared upregulation or downregulation in the sulfidic ecotype compared to nonsulfidic)
- Run linear regression and characterize enhancer-gene pairs that fall within 1 standardized residual of the resulting regression line (i.e., highly correlated enhancer-gene pairs)
- Generate enhancer-gene correlation dot plot for Figure 5. 
- Identify candidate H2S detox genes putatively targeted by a putative enhancer. 

Link to Rscript: [5b_2_putEnhAnalysis_02.02.24.R](https://github.com/blairperry/pmex_csrnaseq/blob/main/analysis/5b_2_putEnhAnalysis_02.02.24.R)

---


