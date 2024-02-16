
library(tidyverse)
library(Biostrings)
library(colorspace)
library(ggrepel)
library(patchwork)
library(org.Hs.eg.db)


# Read in peak tables -----------------------------------------------------

all_peaks <- read_csv('analysis/00_jointAnalyses/new_DI_sig_DE_sig_01.15.24.csv')

# Subset to promoter-TSS peaks
prom_peaks <- all_peaks %>% filter(annotation=='Promoter-TSS' & (logfc_peak * logfc_rna) > 0) # 154 promoter peaks with consistent peak/gene directions

# Split by direction
prom_up <- prom_peaks %>% filter(logfc_peak>0)
prom_dw <- prom_peaks %>% filter(logfc_peak<0)


# Read in candidate gene peaks --------------------------------------------

cand_peaks <- read_tsv('analysis/0_DEanalysis_new/csRNA_DEresults_CandGenes_01.26.24.tsv')
prom_up_cand <- prom_up %>% filter(peak_id %in% cand_peaks$peak_id)

# Read in fasta file  -----------------------------------------------------

all_fasta <- readDNAStringSet('analysis/000_peakTFBS/target.fa')

prom_up_fasta <- all_fasta[names(all_fasta) %in% prom_up$peak_id]
prom_dw_fasta <- all_fasta[names(all_fasta) %in% prom_dw$peak_id]

prom_up_cand_fasta <- all_fasta[names(all_fasta) %in% prom_up_cand$peak_id]

# writeXStringSet(prom_up_fasta,filepath = 'analysis/000_peakTFBS/promoterPeak_Upreg_01.17.24.fa',format = 'fasta')
# writeXStringSet(prom_dw_fasta,filepath = 'analysis/000_peakTFBS/promoterPeak_Downreg_01.17.24.fa',format = 'fasta')
# writeXStringSet(prom_up_cand_fasta,filepath = 'analysis/000_peakTFBS/promoterPeak_CandGeneUpreg_01.26.24.fa',format = 'fasta')


# Get background of non-significant promoter-TSS peaks  --------------------------------
# FDR > 0.5 in both csRNA and mRNA

all_genes <- read_tsv('analysis/0_DEanalysis_new/mRNA_DEResults_All_01.15.24.tsv') %>% janitor::clean_names()
nonsig_genes <- all_genes %>% filter(fdr > 0.05 & abs(log_fc) < 0.5)  


nonsig_peaks <- read_tsv('analysis/0_DEanalysis_new/csRNA_DEResults_All_01.15.24.tsv') %>% 
  janitor::clean_names() %>% 
  filter(str_detect(annotation, 'promoter-TSS')) %>% 
  filter(fdr > 0.05 & abs(log_fc) < 0.5) %>% 
  filter(gene_name %in% nonsig_genes$gene_name)

# 7443 peaks meeting criteria

nonsig_fasta <- all_fasta[names(all_fasta) %in% nonsig_peaks$peak_id]
# writeXStringSet(nonsig_fasta,filepath = 'analysis/000_peakTFBS/nonsigPeak_Background_01.25.24.fa',format = 'fasta')


# Read in CIIIDER results -------------------------------------------------

up_tfbs <- read_csv('analysis/000_peakTFBS/upreg_ciiiderRes/Enrichment: Text9061128605065993800_MostSigDeficit.csv') %>% janitor::clean_names() %>% 
  filter(gene_p_value < 0.01 & gene_representation == 'Up') %>% 
  arrange(-significance_score) %>% 
  mutate(target_proportion = no_transcription_factor_search_genes / total_no_search_genes) %>% 
  rowid_to_column() %>% 
  mutate(peak_set='upregulated')

dw_tfbs <- read_csv('analysis/000_peakTFBS/downreg_ciiiiderRes/Enrichment: Text7404221968844569385_MostSigDeficit.csv') %>% janitor::clean_names() %>% 
  filter(gene_p_value < 0.01 & gene_representation == 'Up') %>% 
  arrange(-significance_score)%>% 
  mutate(target_proportion = no_transcription_factor_search_genes / total_no_search_genes) %>% 
  rowid_to_column() %>% 
  mutate(peak_set = 'downregulated')

nrow(up_tfbs)
nrow(dw_tfbs)

sum(up_tfbs$transcription_factor_id %in% dw_tfbs$transcription_factor_id)


nrow(dw_tfbs) + nrow(up_tfbs)

up <- ggplot(up_tfbs,aes(x=target_proportion,y=log2_enrichment,color=significance_score,size=significance_score)) +
  geom_point(alpha=0.8) +
  geom_text_repel(data=up_tfbs %>% filter(rowid<11),aes(label=transcription_factor_name),size=4,nudge_x = .2,nudge_y = .5,force = 6) +
  labs(x='Percent of peaks with TFBS',y='log2(Enrichment Ratio)',color='-log10(p-value)',size='-log10(p-value)') +
  scale_color_continuous_sequential(palette = 'Viridis') +
  scale_size_binned() +
  theme_linedraw() + theme(panel.grid = element_blank())

dw <- ggplot(dw_tfbs,aes(x=target_proportion,y=log2_enrichment,color=significance_score,size=significance_score)) +
  geom_point(alpha=0.8) +
  geom_text_repel(data=dw_tfbs %>% filter(rowid<11),aes(label=transcription_factor_name),size=4,nudge_x = .2,nudge_y = .5,force = 3) +
  labs(x='Percent of peaks with TFBS',y='log2(Enrichment Ratio)',color='-log10(p-value)',size='-log10(p-value)') +
  scale_color_continuous_sequential(palette = 'Viridis') +
  scale_size_binned() +
  theme_linedraw() + theme(panel.grid = element_blank())

up + dw 

# Plot with top 20 TFBS in each subset ------------------------------------

up_tfbs_top20 <- up_tfbs %>% top_n(n = 20,wt = significance_score) %>% mutate(direction='Parallel Upregulated') 
dw_tfbs_top20 <- dw_tfbs %>% top_n(n = 20,wt = significance_score) %>% mutate(direction='Parallel Downregulated')

all_tfbs_top20 <- up_tfbs_top20 %>% bind_rows(dw_tfbs_top20) %>% mutate(direction = factor(direction,levels=c('Parallel Upregulated','Parallel Downregulated')))

up_top20 <- ggplot(up_tfbs_top20,aes(y=significance_score,x=log2_enrichment,size=target_proportion)) +
  geom_point(alpha=0.8) +
  geom_text_repel(aes(label=transcription_factor_name),size=3) +
  labs(y='-log10(p-value)',x='log2(Enrichment Ratio)',size='Percent of peaks with TFBS') +
  # scale_color_continuous_sequential(palette = 'Viridis') +
  scale_size_continuous(breaks = c(0,.25,.5,.75,1),range = c(2,6)) +
  theme_linedraw() + theme(panel.grid = element_blank())
up_top20

dw_top20 <- ggplot(dw_tfbs_top20,aes(y=significance_score,x=log2_enrichment,size=target_proportion)) +
  geom_point(alpha=0.8) +
  geom_text_repel(aes(label=transcription_factor_name),size=3) +
  labs(y='-log10(p-value)',x='log2(Enrichment Ratio)',size='Percent of peaks with TFBS') +
  # scale_color_continuous_sequential(palette = 'Viridis') +
  scale_size_continuous(breaks = c(0,.25,.5,.75,1),range = c(2,6)) +
  theme_linedraw() + theme(panel.grid = element_blank())

dw_top20

both_top20 <- ggplot(all_tfbs_top20,aes(y=significance_score,x=log2_enrichment,size=target_proportion,color=target_proportion)) +
  geom_point(alpha=0.8) +
  labs(y='-log10(p-value)',x='log2(Enrichment Ratio)',size='Percent of peaks with TFBS',color='Percent of peaks with TFBS') +
  scale_color_binned_sequential(palette = 'Viridis',rev = F,breaks = seq(0,1,by=0.25),limits=c(0,1),end = 1) +
  geom_point(alpha=0.8,pch=21,fill='NA',color='black') +
  geom_text_repel(aes(label=transcription_factor_name),size=3,color='black') +
  facet_wrap(~direction,ncol=2,scales = 'free') +
  scale_size_continuous(breaks = seq(.25,1,by=0.25),range = c(1,5),limits=c(0,1)) +
  # scale_x_continuous(limits = c(0,10)) +
  theme_linedraw() + theme(panel.grid = element_blank(),strip.background = element_rect(fill='NA',color = 'NA'),strip.text = element_text(color='black'))




# Read in and add TF class info -------------------------------------------

tf_info <- read_tsv('analysis/3_tfbs_analysis/tfbs_info/JASPAR2020_CORE_FamilyClass_cleaned.txt') 

up_tfbs.merge <- up_tfbs %>% 
  left_join(tf_info,by=c('transcription_factor_id'='motif_id')) %>% 
  mutate(direction = 'Upregulated Promoter-TSS Peaks')

dw_tfbs.merge <- dw_tfbs %>% 
  left_join(tf_info,by=c('transcription_factor_id'='motif_id')) %>% 
  mutate(direction = 'Downregulated Promoter-TSS Peaks')

all_tfbs.merge <- up_tfbs.merge %>% bind_rows(dw_tfbs.merge) %>% 
  group_by(direction,class) %>% mutate(max_prop = max(target_proportion)) %>% 
  mutate(label = ifelse(target_proportion == max_prop,T,F))


ggplot(up_tfbs.merge,aes(x=target_proportion,y=log2_enrichment,color=class,size=significance_score)) +
  geom_point(alpha=0.8) +
  # geom_text_repel(data=up_tfbs %>% filter(rowid<11),aes(label=transcription_factor_name),size=4,nudge_x = .2,nudge_y = .5,force = 6) +
  labs(x='Percent of peaks with TFBS',y='log2(Enrichment Ratio)',size='-log10(p-value)') +
  scale_size_binned() +
  theme_linedraw() + theme(panel.grid = element_blank())

ggplot(all_tfbs.merge,aes(x=target_proportion,y=class,fill=class,size=significance_score)) +
  geom_point(alpha=0.7,pch=21,show.legend = T) +
  geom_text_repel(data=all_tfbs.merge %>% filter(label==T),aes(label=transcription_factor_name,color=class),size=3,force = 6,show.legend = F) +
  labs(x='Percent of peaks with TFBS',y='TF Class',size='-log10(p-value)') +
  facet_wrap(~direction,ncol = 2) +
  scale_size_continuous(range = c(2,6)) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(0,1),expand = c(0,0)) +
  theme_linedraw() 

ggplot(all_tfbs.merge,aes(x=target_proportion,y=log2_enrichment,fill=class,size=significance_score)) +
  geom_point(alpha=0.7,pch=21,show.legend = T) +
  geom_text_repel(data=all_tfbs.merge %>% filter(rowid<11),aes(label=transcription_factor_name),size=4,nudge_x = .2,nudge_y = .5,force = 6) +
  labs(x='Percent of peaks with TFBS',y='TF Class',size='-log10(p-value)') +
  facet_wrap(~direction,ncol = 2) +
  scale_size_binned() +
  # scale_y_discrete(limits = rev) +
  scale_x_continuous(limits = c(0,1),expand = c(0,0)) +
  theme_linedraw() 








# Identify differentially expressed TFs in mRNAseq dataset ------------------------

# Read in all TFs in Jaspar database and parse gene symbols
all_tf_id <- read_tsv('analysis/000_peakTFBS/jaspar_TFList_01.18.24.txt',col_names = c('motif_id','symbol')) %>% 
  separate_rows(symbol,sep = '[\\::]+') %>% 
  mutate(simple_id=ifelse(
    str_detect(symbol,'\\('), 
    str_split_fixed(symbol,'[(]',2)[,1] %>% str_to_upper(),
    str_to_upper(symbol)
  ))

# Convert symbol to uniprot ID
all_tf.uniprot <- AnnotationDbi::select(org.Hs.eg.db, 
                                   keys = all_tf_id$simple_id,
                                   columns = c("UNIPROT"),
                                   keytype = "SYMBOL") %>% 
  filter(!is.na(UNIPROT)) # remove small number of TFs without a human Uniprot ID

# Add uniprot ID to TF table
all_tf_id.uniprot <- all_tf_id %>% 
  left_join(all_tf.uniprot,by=c('simple_id'='SYMBOL')) # NOTE that some symbols have multiple uniprot IDs, so this table is longer than above

# Read in DE genes and filter to TFs
sig_de <- read_tsv('analysis/0_DEanalysis_new/mRNA_DEResults_Sig_01.15.24.tsv') %>% 
  filter(!is.na(human_id)) %>% 
  mutate(uniprot_id = str_split_fixed(human_id,'[|]',3)[,2]) %>% 
  filter(uniprot_id %in% all_tf_id.uniprot$UNIPROT)

# 103 DE TF genes between ecotypes

sig_de %>% group_by(logFC < 0) %>% tally()
# 47 downregulated, 56 upregulated

# Associate sig TF genes with motif IDs
sig_de.motifs <- sig_de %>% 
  left_join(all_tf_id.uniprot,by=c('uniprot_id'='UNIPROT'))

# Filter to DE TF genes with sig TFs
all_sig_tfs <- up_tfbs %>% bind_rows(dw_tfbs)

sig_de.sig_motifs <- sig_de.motifs %>% 
  filter(motif_id %in% all_sig_tfs$transcription_factor_id) %>% 
  unique() %>% 
  left_join(all_sig_tfs,by=c('motif_id'='transcription_factor_id')) %>% 
  mutate(tf_rna_direction = ifelse(logFC < 0, 'downregulated','upregulated'))

sig_de.sig_motifs %>% group_by(peak_set,tf_rna_direction) %>% tally()

sig_de.sig_motifs %>% dplyr::select(transcription_factor_name) %>% unique() # 26 

same_dir <- sig_de.sig_motifs %>% filter(peak_set == tf_rna_direction)
same_dir %>% filter(tf_rna_direction=='upregulated') %>% 
  dplyr::select(gene.name,human_id,transcription_factor_name) %>% 
  unique() # 9 upregulated TFs that have motifs enriched in the upreg peak/gene set

same_dir %>% filter(tf_rna_direction=='downregulated') %>% 
  dplyr::select(gene.name,human_id,transcription_factor_name) %>% 
  unique() # 3 downregulated TFs have motifs enriched in downreg set (ASCL1 is in there twice)

ggplot(sig_de.sig_motifs %>% filter(tf_rna_direction=='upregulated'),aes(x=logFC,y=log2_enrichment,color=peak_set,size=significance_score)) +
  geom_point() +  
  theme_linedraw()

ggplot(sig_de.sig_motifs %>% filter(tf_rna_direction=='downregulated'),aes(x=logFC,y=log2_enrichment,color=peak_set,size=significance_score)) +
  geom_point() +  
  theme_linedraw()



# Building Supplementary Table S9  ----------------------------------------
# All TF enrichment results

supp_table_s9 <- all_tfbs.merge %>% 
  arrange(rev(direction),transcription_factor_name) %>% 
  dplyr::select(
    `Peak Group`= direction,
    `TF Motif ID`=transcription_factor_id,
    `TF Name` = transcription_factor_name,
    `TF Family` = family,
    `TF Class` = class,
    `Deficit` = deficit,
    `Number of Foreground Peaks with Motif` = no_transcription_factor_search_genes,
    `Number of Background Peaks with Motif` = no_transcription_factor_background_genes,
    `Proportion of Target Sites` = target_proportion,
    `p-value` = gene_p_value,
    `Log2 Enrichment` = log2_enrichment,
    `Significance Score` = significance_score
  ) 

# write_tsv(supp_table_s9,'analysis/000_peakTFBS/suppTableS9_TFenrich_02.02.24.tsv')


# Building Supp Table S10 -------------------------------------------------
# All enriched TF binding sites with same direction DE

names(same_dir)
supp_table_s10 <- same_dir %>% 
  dplyr::select(
    `TF Motif ID` = motif_id,
    `TF Name` = transcription_factor_name,
    `Log2 Enrichment Ratio` = log2_enrichment,
    `Enrichment p-value` = gene_p_value,
    `Enhancer-Gene Direction` = peak_set,
    `TF Encoding Gene ID` = gene.name,
    `Human Homologue`=human_id,
    `TF Encoding Gene Log2 Fold Change` = logFC,
    `TF Encoding Gene p-value` = PValue,
    `TF Encoding Gene FDR` = FDR
  )


write_tsv(supp_table_s10,'analysis/000_peakTFBS/suppTableS10_TFGeneExp_02.02.24.tsv')


#
##