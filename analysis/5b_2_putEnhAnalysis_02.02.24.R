
library(tidyverse)
library(cols4all)

# Read in csRNAseq and RNAseq DE results ----------------------------------

peak.all <- read_tsv('analysis/0_DEanalysis_new/csRNA_DEResults_All_01.15.24.tsv') %>% janitor::clean_names()
rna.all <- read_tsv('analysis/0_DEanalysis_new/mRNA_DEResults_All_01.15.24.tsv')  %>% janitor::clean_names()
  
peak.de <- peak.all %>% filter(fdr < 0.05)
nrow(peak.de)

rna.de <- rna.all %>% filter(fdr < 0.05)
rna.de %>% nrow()
rna.all.notDe <- rna.all %>% filter(fdr > 0.05)

# All putative enhancers
put_enh <- read_tsv('analysis/1_peakExploration/putativeEnhancerPeaks_pmex_01.12.24.bed',
                    col_names = c('chr','start','end','id'))

put_enh %>% filter(id %in% peak.de$peak_id) %>% nrow()
# 448

# Assess genes within distance threshold ------------------

colnames <- c('chr_peak','startSlop_peak','endSlop_peak','id_peak','ignore1','start_tss','end_tss','desc')

genes.500kb <- read_tsv('analysis/1_PeakExploration/putativeEnhancerPeaks_pmex_01.12.24.500kbSlop.tssIntersect.txt',col_names = colnames) %>% select(-contains('ignore')) %>% mutate(gene = str_split_fixed(desc,'[;]',2)[,1] %>% str_remove_all('ID=gene-'))

genes.500kb.de <- genes.500kb %>% filter(id_peak %in% peak.de$peak_id)

# Calculate percent of DE genes within X distance of ANY enhancers
(rna.de %>% filter(gene_name %in% genes.500kb$gene) %>% nrow()) / nrow(rna.de) # ~90%

# Calculate percent of DE genes within X distance of DE enhancers
(rna.de %>% filter(gene_name %in% genes.500kb.de$gene) %>% nrow()) / nrow(rna.de) # ~19%


# Merge DE genes and DE enhancer tables
rna.peak.de.500kb <- rna.de %>% filter(gene_name %in% genes.500kb.de$gene) %>% 
  left_join(genes.500kb.de,by=c('gene_name'='gene')) %>% 
  mutate(center_peak = str_split_fixed(id_peak,'[-]',4)[,3] %>% as.numeric()) %>% 
  mutate(peak_tss_absDist = abs(start_tss - center_peak)) %>% 
  left_join(peak.de,by=c('id_peak'='peak_id'))

rna.peak.notDe.500kb <- rna.all.notDe %>% filter(gene_name %in% genes.500kb.de$gene) %>% 
  left_join(genes.500kb.de,by=c('gene_name'='gene')) %>% 
  mutate(center_peak = str_split_fixed(id_peak,'[-]',4)[,3] %>% as.numeric()) %>% 
  mutate(peak_tss_absDist = abs(start_tss - center_peak))

rna.peak.de.500kb %>% group_by(gene_id) %>% tally() %>% summarise(mean= mean(n), median = median(n)) # mean 1.38, median 1 DE enhancers per DE gene

ggplot(rna.peak.de.500kb,aes(x=log_fc.x,y=log_fc.y)) +
  geom_point() +
  theme_linedraw()

rna.peak.de.500kb %>% 
  mutate(rna_direction = ifelse(log_fc.x > 0,'Up','Down')) %>% 
  mutate(peak_direction = ifelse(log_fc.y > 0,'Up','Down')) %>% 
  group_by(rna_direction,peak_direction) %>% 
  tally()

# # Groups:   rna_direction [2]
# rna_direction peak_direction     n
# <chr>         <chr>          <int>
# 1 Down          Down             238
# 2 Down          Up               161
# 3 Up            Down             164
# 4 Up            Up               177

# 415 peak/RNA DE in same direction
# 325 peak/RNA DE in opposite direction


# Positively correlated - putative positive regulators ------------------------------

rna.peak.de.500kb.posCorr <- rna.peak.de.500kb %>% 
  dplyr::rename(rna_logfc = log_fc.x, csrna_logfc = log_fc.y) %>% 
  mutate(rna_direction = ifelse(rna_logfc > 0,'Up','Down')) %>% 
  mutate(peak_direction = ifelse(csrna_logfc > 0,'Up','Down')) %>% 
  filter(rna_direction == peak_direction)

ggplot(rna.peak.de.500kb.posCorr,aes(x=rna_logfc,y=csrna_logfc)) +
  geom_point() +
  theme_linedraw()

model.pos <- lm(csrna_logfc ~ rna_logfc, data=rna.peak.de.500kb.posCorr)
model.sum.pos <- summary(model.pos)

model.sum.pos
# y = 1.307203x + 00.07021
# R2 = 0.6179
# p < 2.2e-16

#calculate the standardized residuals
standard_res.pos <- rstandard(model.pos)

#view the standardized residuals
standard_res.pos

rna.peak.de.500kb.posCorr.stdRes <- rna.peak.de.500kb.posCorr %>% 
  mutate(stdRes = standard_res.pos) %>% 
  mutate(stdRes_threshold = ifelse(abs(stdRes) > 1,'out','in')) %>% 
  mutate(quadrant = case_when(
    stdRes_threshold=='in' & rna_logfc>0~'positive',
    stdRes_threshold=='in' & rna_logfc<0~'negative',
    stdRes_threshold=='out' ~ 'out'
  ))

rna.peak.de.500kb.posCorr.stdRes %>% group_by(stdRes_threshold) %>% tally()

enh_dotplot <- ggplot(rna.peak.de.500kb.posCorr.stdRes,aes(x=rna_logfc,y=csrna_logfc,color=quadrant)) +
  geom_point(alpha=0.7,show.legend = F) +
  geom_abline(slope = model.sum.pos$coefficients[2,1],
              intercept = model.sum.pos$coefficients[1,1],
              lty=2) +
  labs(x='log2FoldChange - Gene RNAseq',y='log2FoldChange - Peak csRNAseq', color= 'Standard Residual Threshold = 1') +
  scale_color_manual(values = c('positive'='firebrick','negative'='darkblue','out'='grey70')) +
  coord_cartesian(xlim=c(-6,6),ylim=c(-6,6)) +
  geom_vline(xintercept = 0,lwd=0.25) +
  geom_hline(yintercept = 0,lwd=0.25) +
  theme_linedraw() + theme(panel.grid = element_blank())

enh_dotplot


# Filter to peaks within standard residual threshold
rna.peak.de.500kb.allCorr.stdRes.inThresh <- rna.peak.de.500kb.posCorr.stdRes %>% filter(stdRes_threshold!='out') %>% 
  mutate(peak_id = paste(chr,start,end,sep='_'))

rna.peak.de.500kb.allCorr.stdRes.inThresh %>% select(peak_id) %>% unique() %>% nrow() # 213 PRRs
rna.peak.de.500kb.allCorr.stdRes.inThresh %>% select(gene_id) %>% unique() %>% nrow() # 286 genes

rna.peak.de.500kb.allCorr.stdRes.inThresh %>% 
  mutate(rna_dir = ifelse(rna_logfc > 0, 'positive','negative'),
         csrna_dir = ifelse(csrna_logfc > 0, 'positive','negative')) %>% 
  select(peak_id,rna_dir,csrna_dir) %>% 
  unique() %>% 
  group_by(rna_dir,csrna_dir) %>% 
  tally()
# 86 unique upregulated PRRs associated with upregulated genes
# 127 unique downregulated PRRs associated with downregulated genes


rna.peak.de.500kb.allCorr.stdRes.inThresh %>% # histogram of genes per peak
  group_by(peak_id) %>% 
  tally() %>% 
  ggplot(aes(x=n)) +
  geom_histogram() + theme_linedraw()

rna.peak.de.500kb.allCorr.stdRes.inThresh %>%
  group_by(peak_id) %>% 
  tally() %>% 
  summarise(mean=mean(n),median=median(n))

rna.peak.de.500kb.allCorr.stdRes.inThresh


# Build Supp Table S11 - all putative enhancer-gene pairs -----------------
names(rna.peak.de.500kb.allCorr.stdRes.inThresh)
supp_table_s11 <- rna.peak.de.500kb.allCorr.stdRes.inThresh %>% 
  select(`csRNAseq Peak ID` = id_peak,
         Chromosome = chr,
         `csRNAseq Peak Start` = start,
         `csRNAseq Peak End` = end,
         `csRNAseq Strand` = strand,
         `csRNA-seq log2 Fold Change` = csrna_logfc,
         `csRNA-seq FDR` = fdr.y,
         `Inferred Target Gene`=gene_name.x,
         `Human Homologue` = human_id,
         `Human Gene Description` = protein_annotations,
         `Target Gene Log2 Fold Change`=rna_logfc,
         `Target Gene FDR` = fdr.x,
         `Peak-to-Gene Absolute Distance (bp)`=peak_tss_absDist)

# write_tsv(supp_table_s11,'analysis/2_csRNA_RNAseq_Integration/SuppTableS11_EnhGenePairs_02.05.24.tsv')


# write_tsv(rna.peak.de.500kb.allCorr.stdRes.inThresh,'analysis/2_csRNA_RNAseq_Integration/pmex_PERs_Corr_AllInfoUPDATE_02.02.24.tsv')

rna.peak.de.500kb.allCorr.stdRes.inThresh %>% 
  mutate(length = abs(start - end)) %>% 
  select(peak_id,length) %>% unique() %>% 
  summarise(meanLen = mean(length)) 
# Mean length 167

ggplot(rna.peak.de.500kb.allCorr.stdRes.inThresh,aes(x=peak_tss_absDist)) + 
  geom_density(fill='grey',alpha=0.5) +
  theme_linedraw()

# Make bed of correlated peaks
rna.peak.de.500kb.allCorr.stdRes.inThresh.bed <- rna.peak.de.500kb.allCorr.stdRes.inThresh %>% 
  mutate(length = abs(start - end)) %>% 
  select(id_peak,chr,start,end,length) %>% 
  mutate(peak_id = paste(id_peak,'_len',length,sep = '')) %>% 
  select(chr,start,end,peak_id) %>% 
  unique() %>% 
  arrange(chr,start)
  
# write_tsv(rna.peak.de.500kb.allCorr.stdRes.inThresh.bed,'analysis/2_csRNA_RNAseq_Integration/pmex_PERs_Corr_ActualSize_01.15.24.bed',col_names = F)

# Output list of all genes with correlated peaks
rna.peak.de.500kb.allCorr.stdRes.inThresh.genes <- rna.peak.de.500kb.allCorr.stdRes.inThresh %>% 
  select(gene_id,gene_name.x,human_id,rna_logfc,rna_fdr=fdr.x) %>% 
  unique()

corr.genes.up <- rna.peak.de.500kb.allCorr.stdRes.inThresh.genes %>% filter(rna_logfc>0) %>% # 110
  select(gene_id,gene_name.x,human_id,rna_logfc,rna_fdr) %>% 
corr.genes.dw <- rna.peak.de.500kb.allCorr.stdRes.inThresh.genes %>% filter(rna_logfc<0) %>% # 176
  select(gene_id,gene_name.x,human_id,rna_logfc,rna_fdr)

# write_tsv(corr.genes.up,'analysis/2_csRNA_RNAseq_Integration/pmex_corrGenes_up_UPDATE_02.02.24.tsv')
# write_tsv(corr.genes.dw,'analysis/2_csRNA_RNAseq_Integration/pmex_corrGenes_dw_UPDATE_02.02.24.tsv')


# Test overlap with candidate gene set ------------------------------------

# Merge significant peaks with H2S-related candidate list. The Sulfide Detox/Response Gene Set is relevant for this study. 
# The NuclearRef and Broughton Gene Sets are reference sets.

# Read in candidate list
cand_list <- read.csv('data_fromKerryAnalysis/reference/OXPHOS_Reference_and_Detox_Gene_IDS.csv', header = 1)

# Subset to only include Sulfide Detox/Response Gene Set
subset_cand <- subset(cand_list, Gene.Set == "Sulfide Detox/Response")

# Change name of sulfide:quinone oxidoreductase from "sqor" to "sqrdl" in Gene.name.from.accession column
subset_cand[subset_cand == "sqor"] <- "sqrdl"

# Merge H2S candidate list with up and downregulated genes with correlated enhancer
merged_up <- corr.genes.up %>% left_join(subset_cand,by=c('gene_name.x'='Gene.name.from.accession')) %>%
  filter(!is.na(Gene.Set))

merged_up_peak <- rna.peak.de.500kb.allCorr.stdRes.inThresh %>% filter(gene_name.x %in% subset_cand$Gene.name.from.accession)

# write_tsv(merged_up_peak,'analysis/2_csRNA_RNAseq_Integration/pmex_corrGenes_CandGenes_UPDATE_02.02.24.tsv')

# Determine if the number of differentially expressed peaks in the Sulfide Detox/Response Gene Set is significant using a Fisher's Exact Test.

#2 x 2 contingency table (Set = Sulfide Detox/Response Gene Set):

### Upreg genes targeted by DE enhancers

up_target_cand <- nrow(merged_up)
up_target_notcand <- nrow(corr.genes.up) - nrow(merged_up)

notTarget_cand <- nrow(subset_cand) - nrow(merged_up)
notTarget_notcand <- rna.all %>% 
  filter(!(gene_name %in% subset_cand$Gene.name.from.accession)) %>% 
  filter(!(gene_name %in% corr.genes.up$gene_name.x)) %>% 
  nrow()

# Make the 2 x 2 contingency table as a matrix
table_up <- matrix(c(up_target_cand, up_target_notcand, notTarget_cand, notTarget_notcand),
                nrow = 2,
                dimnames = list(c('in_set', 'not_in_set'),
                                c('target_up', 'not_target_up')))
# Fisher's Exact Test for count data
fisher.test(table_up, alternative = "two.sided") # Sig, p=0.2307

# NOT ENRICHED - two genes (slc26a5, slc26a11)

