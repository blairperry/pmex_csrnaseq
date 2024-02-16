
library(tidyverse)
library(org.Hs.eg.db)
library(ggrepel)
library(colorspace)

# Read in all peaks, enhancers, and gene-associated enhancers -------------------------------------------------------
peak.all <- read_tsv('analysis/0_DEanalysis_new/csRNA_DEResults_All_01.15.24.tsv',comment = '#') %>% janitor::clean_names()
enh.all <- read_tsv('analysis/1_peakExploration/putativeEnhancerPeaks_pmex_01.12.24.bed',col_names = c('chr','start','end','peak_id'))
corr_enh.all <- read_tsv('analysis/2_csRNA_RNAseq_Integration/pmex_PERs_Corr_AllInfo_01.15.24.tsv')

cand_enh <- read_tsv('analysis/2_csRNA_RNAseq_Integration/pmex_corrGenes_CandGenes_01.26.24.tsv')

# Get DE gene-associated enhancer peaks ------------------------------------------------------------

corr_enh.all_sameDir <- corr_enh.all %>% filter((rna_logfc * csrna_logfc) > 0)

corr_enh.de <- peak.all %>% 
  filter(peak_id %in% corr_enh.all_sameDir$id_peak) %>% 
  filter(fdr < 0.05)
nrow(corr_enh.de)

corr_enh.de.up <- corr_enh.de %>% filter(log_fc > 0)
corr_enh.de.dw <- corr_enh.de %>% filter(log_fc < 0)

corr_enh.de.up.bed <- enh.all %>% filter(peak_id %in% corr_enh.de.up$peak_id)%>% 
  arrange(chr,start,end)
corr_enh.de.dw.bed <- enh.all %>% filter(peak_id %in% corr_enh.de.dw$peak_id)%>% 
  arrange(chr,start,end)

cand_enh.de.up.bed <- corr_enh.all %>% filter(peak_id %in% cand_enh$peak_id)%>% 
  dplyr::select(chr,start,end,peak_id) %>% 
  unique() %>% 
  arrange(chr,start,end)


# Get 1000 non-DE putative enhancers for background -----------------------
set.seed(1)

bg_enh <- peak.all %>% 
  filter(peak_id %in% enh.all$peak_id) %>% 
  filter(fdr > 0.05) %>% 
  sample_n(1000)

bg_enh.bed <- enh.all %>% filter(peak_id %in% bg_enh$peak_id) %>% 
  arrange(chr,start,end)


# Write bedfiles  ---------------------------------------------------------

# write_tsv(corr_enh.de.up.bed,'analysis/3_tfbs_analysis/corr_enhancers_upreg_01.19.24.bed',col_names = F)
# write_tsv(corr_enh.de.dw.bed,'analysis/3_tfbs_analysis/corr_enhancers_dwreg_01.19.24.bed',col_names = F)
# write_tsv(bg_enh.bed,'analysis/3_tfbs_analysis/background_enhancers_01.19.24.bed',col_names = F)

# write_tsv(cand_enh.de.up.bed,'analysis/3_tfbs_analysis/candGene_enhancers_upreg_01.26.24.bed',col_names = F)




# Read in CIIIDER results  ------------------------------------------------

upreg_tfs <- read_csv('analysis/3_tfbs_analysis/ciiider_results/upreg/Enrichment: Text7881396120718037618_MostSigDeficit.csv') %>% janitor::clean_names()
dwreg_tfs <- read_csv('analysis/3_tfbs_analysis/ciiider_results/downreg/Enrichment: Text2007677463292186147_MostSigDeficit.csv') %>% janitor::clean_names()

upreg_tfs.sigEnrich <- upreg_tfs %>% filter(gene_p_value < 0.01 & gene_representation=='Up') %>% 
  arrange(-significance_score)%>% 
  mutate(target_proportion = no_transcription_factor_search_genes / total_no_search_genes) %>% 
  rowid_to_column() %>% 
  mutate(peak_set = 'upregulated',
         peak_group = 'Upregualted Enhancer-Gene Pairs')

dwreg_tfs.sigEnrich <- dwreg_tfs %>% filter(gene_p_value < 0.01 & gene_representation=='Up') %>% 
  arrange(-significance_score)%>% 
  mutate(target_proportion = no_transcription_factor_search_genes / total_no_search_genes) %>% 
  rowid_to_column() %>% 
  mutate(peak_set = 'downregulated',
         peak_group = 'Downregulated Enhancer-Gene Pairs')

sum(upreg_tfs.sigEnrich$transcription_factor_name %in% dwreg_tfs.sigEnrich$transcription_factor_name)

ggplot(upreg_tfs,aes(x=average_log2_proportion_bound,y=log2_enrichment,color=abs(significance_score))) +
  geom_point() + 
  theme_linedraw()


# Build supp table S12 - all enhancer TF results --------------------------
all_tfbs <- upreg_tfs.sigEnrich %>% 
  bind_rows(dwreg_tfs.sigEnrich)

tf_info <- read_tsv('analysis/3_tfbs_analysis/tfbs_info/JASPAR2020_CORE_FamilyClass_cleaned.txt') 

all_tfbs.info <- all_tfbs %>% 
  left_join(tf_info,by=c('transcription_factor_id'='motif_id'))


supp_table_s12 <- all_tfbs.info %>% 
  arrange(rev(peak_group),transcription_factor_name) %>% 
  dplyr::select(
    `Peak Group`= peak_group,
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

write_tsv(supp_table_s12,'analysis/3_tfbs_analysis/suppTableS12_EnhTFenrich_02.05.24.tsv')

##



# Plots of top 20 enriched TFBS -------------------------------------------



up_tfbs_top20 <- upreg_tfs.sigEnrich %>% top_n(n = 20,wt = significance_score) %>% mutate(direction='Parallel Upregulated') 
dw_tfbs_top20 <- dwreg_tfs.sigEnrich %>% top_n(n = 20,wt = significance_score) %>% mutate(direction='Parallel Downregulated')

all_tfbs_top20 <- up_tfbs_top20 %>% bind_rows(dw_tfbs_top20) %>% mutate(direction = factor(direction,levels=c('Parallel Upregulated','Parallel Downregulated')))


both_top20 <- ggplot(all_tfbs_top20,aes(y=significance_score,x=log2_enrichment,size=target_proportion,color=target_proportion)) +
  geom_point(alpha=0.8) +
  labs(y='-log10(p-value)',x='log2(Enrichment Ratio)',size='Percent of peaks with TFBS',color='Percent of peaks with TFBS') +
  scale_color_binned_sequential(palette = 'Viridis',rev = F,breaks = seq(0,1,by=0.25),limits=c(0,1),end = 1) +
  geom_point(alpha=0.8,pch=21,fill='NA',color='black') +
  geom_text_repel(aes(label=transcription_factor_name),size=3,color='black') +
  facet_wrap(~direction,ncol=2,scales = 'free') +
  scale_size_continuous(breaks = seq(.25,1,by=0.25),range = c(1,5),limits=c(0,1)) +
  # scale_x_continuous(limits = c(0,10)) +
  theme_linedraw() + theme(panel.grid = element_blank())

# Only works if enh_dotplot from DEIntegration... Rscript is in environment
enh_dotplot + both_top20 + plot_layout(widths = c(1,2))

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

# Associate sig TF genes with motif IDs
sig_de.motifs <- sig_de %>% 
  left_join(all_tf_id.uniprot,by=c('uniprot_id'='UNIPROT'))

# Filter to DE TF genes with sig TFs
all_sig_tfs <- upreg_tfs.sigEnrich %>% bind_rows(dwreg_tfs.sigEnrich)

sig_de.sig_motifs <- sig_de.motifs %>% 
  filter(motif_id %in% all_sig_tfs$transcription_factor_id) %>% 
  unique() %>% 
  left_join(all_sig_tfs,by=c('motif_id'='transcription_factor_id')) %>% 
  mutate(tf_rna_direction = ifelse(logFC < 0, 'downregulated','upregulated'))

sig_de.sig_motifs %>% group_by(peak_set,tf_rna_direction) %>% tally()


same_dir <- sig_de.sig_motifs %>% filter(peak_set == tf_rna_direction)

same_dir %>% filter(tf_rna_direction=='upregulated') %>% 
  dplyr::select(gene.name,human_id,transcription_factor_name) %>% 
  unique() 

same_dir %>% filter(tf_rna_direction=='downregulated') %>% 
  dplyr::select(gene.name,human_id,transcription_factor_name) %>% 
  unique() 



# Build Supp table S13 - enh TFs with same direction DE  ------------------


names(same_dir)
supp_table_s13 <- same_dir %>% 
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

write_tsv(supp_table_s13,'analysis/3_tfbs_analysis/suppTableS13_EnhTFGeneExp_02.05.24.tsv')

#
##


#




