# Analysis of transcription in H2S-adapted fish (*Poecilia mexicana*)

This repo contains code for all data processing and analysis used in our study of transcription and gene expression in H2S-adapted 
populations of *P. mexicana*.

For questions, please contact blair.perry(at)wsu.edu.

## Data Availability
- Newly-generated csRNA-seq data
	- Raw data available at NCBI BioProject: PRJNA743555
- Previously published mRNA-seq data ([Kelley et al. 2016](http://mbe.oxfordjournals.org/content/early/2016/03/01/molbev.msw020))
	- Raw data available at NCBI BioProject: [PRJNA290391](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA290391/)

## Contents
1. cRNA-seq trimming, mapping, peak calling, and peak annotation
2. mRNA-seq trimming and mapping 
3. Analyses of differential transcript initiation (csRNA) and gene expression (mRNA)
4. Joint analyses of csRNA and mRNA data
	1. Characterization of transcript initiation annotations
	2. Identification of putative enhancer regions
5. TF binding site enrichment analyses


---

### 1. cRNA-seq QC, mapping, peak calling and annotation

These analyses were run on WSU's HPC ([Kamiak](https://hpc.wsu.edu/)). For generalizability, simplified commands are presented here 
rather than the specific SLURM scripts used to run these commands on Kamiak. 

Note: The variable `$sample` should be replaced with the full sample name.
#### Trim with homerTools *trim*
```bash
homerTools trim \
  -3 AGATCGGAAGAGCACACGTCT \
  -mis 2 \
  -minMatchLength 4 \
  -min 20 \
  -q 20 \
  $sample.fastq.gz
```

#### Map csRNA-seq with STAR
Reference genome and annotations in GFF format were downloaded from NCBI (Accession: 
[GCF_001443325.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001443325.1/)).

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

#### Make tag directories
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

#### Identify transcription initiation peaks
Peaks are clusters of transcription initiation, called using the "full approach" of the HOMER pipeline. This incorporates mRNA tag 
directories into the csRNA peak calling to help eliminate false positives like exomic contaminants.

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

#### Merge peaks
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

#### Annotate peaks
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

### 2. mRNA-seq trimming and mapping 
The variable $sample should be replaced with the full sample name.
#### Trim with TrimGalore
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


#### Map RNA-seq with HISAT2

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

#### Convert, sort, and merge SAM/BAM files

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

#### Build gene count matrices

Note: modified prepDE.py script (prepDE-2019-10-01_edited.py) is provided in this repository. 

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

