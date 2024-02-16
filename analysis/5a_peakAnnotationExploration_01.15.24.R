
library(tidyverse)
library(ggforce)

# Intersect Differentially Initiated Transcripts and Differentially Expressed Genes
# Differentially Initiated (DI) peaks from csRNA-seq
# Read in DI peaks (csRNA-seq) analyzed using edgeR.

DI <- read_tsv("analysis/0_DEanalysis_new/csRNA_DEResults_All_01.15.24.tsv") %>% 
  select(-contains('normalization')) %>% 
  janitor::clean_names()
names(DI)


# Split Annotation column, this is needed to later summarize the annotations of the peaks
DI_split <- DI %>% 
  mutate(anno_refseq_accession=str_split_fixed(annotation,' ',2)[,2]) %>% 
  mutate(annotation=str_split_fixed(annotation,' ',2)[,1])

# Capitalize the first letter of all annotations
str_sub(DI_split$annotation, 1, 1) <- str_sub(DI_split$annotation, 1, 1) %>% str_to_upper()

#Sort by FDR (adjusted p-value)
DI_sorted <- DI_split %>% arrange(fdr)

# Subset to significantly DI peaks only (FDR < 0.05)
DI_sig <- DI_sorted %>% filter(fdr < 0.05)

print(paste0("Number of significantly DI peaks: ", nrow(DI_sig)))

DI_sig %>% filter(log_fc>0) %>% nrow() # 550 upreg
DI_sig %>% filter(log_fc<0) %>% nrow() # 511 downreg
nrow(DI_sig)


# Summarize the annotations of all significantly DI peaks (i.e., what region of the genome the peaks map to).

# Count unnanotated DI peaks (NAs)
print(paste0("Number of unannotated peaks: ", sum(is.na(DI_sig$annotation))))
DI_sig.noNA <- DI_sig %>% filter(!is.na(annotation)) # remove peaks that were not associated with gene

# Sum occurrences by genomic peak location
DI_sig_gen_loc <- DI_sig.noNA %>% 
  group_by(annotation) %>% 
  tally() %>% 
  mutate(perc = n / sum(n))

DI_sig_gen_loc  

# Differentially Expressed (DE) genes from RNA-seq
# Read in DE genes (RNA-seq) analyzed using edgeR.

DE <- read_tsv("analysis/0_DEanalysis_new/mRNA_DEResults_All_01.15.24.tsv") %>% janitor::clean_names()

# Subset to significantly DE genes only (FDR < 0.05)
DE_sig <- DE %>% filter(fdr<0.05)
print(paste0("Number of significantly DE genes: ", nrow(DE_sig)))


# Intersect significantly DI peaks with DE genes
# Note, we start with only significantly DI peaks but all DE genes, regardless of their significance.

all_DI_DE <- merge(x = DI_sig.noNA, y = DE, by.x = "gene_name", by.y = "gene_name", all = TRUE)
names(all_DI_DE)

# Subset columns in an order that makes sense
all_DI_DE_subset <- all_DI_DE %>% dplyr::select(gene_id, peak_id, chr, start, end, strand, peak_score, annotation, anno_refseq_accession, distance_to_tss, entrez_id, human_id, nearest_promoter_id, nearest_unigene, logfc_peak=log_fc.x, pval_peak=p_value.x, fdr_peak=fdr.x, logfc_rna=log_fc.y, pval_rna=p_value.y, fdr_rna=fdr.y, protein_annotations)


# Subset significantly DI peaks to significantly DE genes.

DI_sig_DE_sig <- all_DI_DE_subset %>% filter(fdr_peak < 0.05 & fdr_rna < 0.05)
print(paste0("Number of significantly DI peaks with significant DE: ", nrow(DI_sig_DE_sig)))
# write_csv(DI_sig_DE_sig, "analysis/00_jointAnalyses/new_DI_sig_DE_sig_01.15.24.csv")

# Get percentage of peaks mapped to DE genes, non-DE genes
de_match <- nrow(DI_sig_DE_sig)
nonDE_match <- all_DI_DE_subset %>% filter(!is.na(peak_id))  %>% filter(fdr_rna>0.05) %>% nrow()
de_match / (de_match + nonDE_match)


# Summarize the annotations of significantly DI peaks that have significant DE (i.e., what region of the genome the peaks map to).

# Sum occurrences by genomic peak location
DI_sig_DE_sig_gen_loc <- as.data.frame(table(DI_sig_DE_sig$annotation))

# Change column header
names(DI_sig_DE_sig_gen_loc)[1] <- "annotation"

# Add column of percentages
DI_sig_DE_sig_gen_loc_w_percents <- DI_sig_DE_sig_gen_loc %>% 
  mutate(Perc = `Freq` / sum(`Freq`)) %>% 
  arrange(Perc) %>% 
  mutate(Labels = scales::percent(Perc))
# Add column with data frame name
DI_sig_DE_sig_gen_loc_w_percents$df <- "DE Genes"
DI_sig_DE_sig_gen_loc_w_percents

# G-test (log likelihood ratio) of goodness of fit for significantly DI peaks that have significant DE. Do the annotations fit the distribution of all DI peaks?
library('DescTools')

# Expected **percentages**
exp_tot <- sum(DI_sig_gen_loc$n)
exp_exon <- DI_sig_gen_loc[DI_sig_gen_loc$annotation=='Exon',2] %>% as.double() / exp_tot
exp_intergenic <- DI_sig_gen_loc[DI_sig_gen_loc$annotation=='Intergenic',2] %>% as.double() / exp_tot
exp_intron <- DI_sig_gen_loc[DI_sig_gen_loc$annotation=='Intron',2] %>% as.double() / exp_tot
exp_prom <- DI_sig_gen_loc[DI_sig_gen_loc$annotation=='Promoter-TSS',2] %>% as.double() / exp_tot
exp_tts <- DI_sig_gen_loc[DI_sig_gen_loc$annotation=='TTS',2] %>% as.double() / exp_tot

expected_all_DI_peaks <- c(exp_exon,exp_intergenic,exp_intron,exp_prom,exp_tts)
expected_all_DI_peaks

# Observed **frequencies**
obs_exon <- DI_sig_DE_sig_gen_loc_w_percents[DI_sig_DE_sig_gen_loc_w_percents$annotation=='Exon',2] %>% as.double() 
obs_intergenic <- DI_sig_DE_sig_gen_loc_w_percents[DI_sig_DE_sig_gen_loc_w_percents$annotation=='Intergenic',2] %>% as.double() 
obs_intron <- DI_sig_DE_sig_gen_loc_w_percents[DI_sig_DE_sig_gen_loc_w_percents$annotation=='Intron',2] %>% as.double() 
obs_prom <- DI_sig_DE_sig_gen_loc_w_percents[DI_sig_DE_sig_gen_loc_w_percents$annotation=='Promoter-TSS',2] %>% as.double()
obs_tts <- DI_sig_DE_sig_gen_loc_w_percents[DI_sig_DE_sig_gen_loc_w_percents$annotation=='TTS',2] %>% as.double() 

observed_DI_sig_DE_sig <- c(obs_exon,obs_intergenic,obs_intron,obs_prom,obs_tts)

GTest(x = observed_DI_sig_DE_sig,
      p = expected_all_DI_peaks,
      correct = "none")
# G = 27.796, X-squared df = 4, p-value = 1.372e-05



# Subset significantly DI peaks to non-significantly DE genes.

DI_sig_DE_nonsig <- all_DI_DE_subset %>% filter(fdr_peak < 0.05 & fdr_rna > 0.05)
print(paste0("Number of significantly DI peaks without significant DE: ", nrow(DI_sig_DE_nonsig)))
# write.csv(x = DI_sig_DE_nonsig, file = "new_DI_sig_DE_nonsig.csv")


# Summarize the annotations of significantly DI peaks that have non-significant DE (i.e., what region of the genome the peaks map to).

# Sum occurrences by genomic peak location
DI_sig_DE_nonsig_gen_loc <- as.data.frame(table(DI_sig_DE_nonsig$annotation))
# Change column header
names(DI_sig_DE_nonsig_gen_loc)[1] <- "annotation"
# Add column of percentages
DI_sig_DE_nonsig_gen_loc_w_percents <- DI_sig_DE_nonsig_gen_loc %>% 
  mutate(Perc = `Freq` / sum(`Freq`)) %>% 
  arrange(Perc) %>% 
  mutate(Labels = scales::percent(Perc))
# Add column with data frame name
DI_sig_DE_nonsig_gen_loc_w_percents$df <- "Not DE"
DI_sig_DE_nonsig_gen_loc_w_percents


# G-test (log likelihood ratio) of goodness of fit for significantly DI peaks that have non-significant DE. Do the annotations fit the distribution of all DI peaks?

expected2_all_DI_peaks <- expected_all_DI_peaks

# In the order: Exon, Intergenic, Intron, Promoter-TSS, TTS
# Observed frequencies
obs2_exon <- DI_sig_DE_nonsig_gen_loc_w_percents[DI_sig_DE_nonsig_gen_loc_w_percents$annotation=='Exon',2] %>% as.double() 
obs2_intergenic <- DI_sig_DE_nonsig_gen_loc_w_percents[DI_sig_DE_nonsig_gen_loc_w_percents$annotation=='Intergenic',2] %>% as.double() 
obs2_intron <- DI_sig_DE_nonsig_gen_loc_w_percents[DI_sig_DE_nonsig_gen_loc_w_percents$annotation=='Intron',2] %>% as.double() 
obs2_prom <- DI_sig_DE_nonsig_gen_loc_w_percents[DI_sig_DE_nonsig_gen_loc_w_percents$annotation=='Promoter-TSS',2] %>% as.double()
obs2_tts <- DI_sig_DE_nonsig_gen_loc_w_percents[DI_sig_DE_nonsig_gen_loc_w_percents$annotation=='TTS',2] %>% as.double() 

observed2_DI_sig_DE_nonsig <- c(obs2_exon,obs2_intergenic,obs2_intron,obs2_prom,obs2_tts)


observed_DI_sig_DE_nonsig <- c(87, 81, 173, 154, 19)
GTest(x = observed2_DI_sig_DE_nonsig,
      p = expected2_all_DI_peaks,
      correct = "none")
# G = 8.027, X-squared df = 4, p-value = 0.0906 ## NOT SIG


# Generate a stacked bar chart of genomic locations (using a colorblind friendly palette).

# Combine data frames
gen_loc_w_percents <- rbind(DI_sig_DE_sig_gen_loc_w_percents, DI_sig_DE_nonsig_gen_loc_w_percents)
gen_loc_w_percents$df <- factor(gen_loc_w_percents$df,levels=c('Not DE','DE Genes'))

# Color palette
colorBlind   <- c("#E69F00", "#56B4E9", "#D55E00", "#009E73", "#CC79A7", "#F0E442", "#F0E442")

# Plot
stacked_bar <- ggplot(data = gen_loc_w_percents, aes(fill = annotation, x = df, y = Perc)) +
  # Stacked bar chart, use position = "stack" for frequency and position = "fill" for percentages
  geom_bar(position = "fill", stat = "identity", width = 0.8,show.legend = F) +
  #scale_fill_viridis(discrete = TRUE, option = "cividis") +
  scale_fill_manual(values = colorBlind) +
  xlab("Intersection") +
  ylab("Frequency") +
  scale_y_reverse() +
  theme(aspect.ratio = .1) +
  theme_classic(base_size = 14) +
  # Flip x and y axes
  coord_flip()
stacked_bar

# Subset significantly DI peaks with significantly DE genes by the directionality of their log2-fold change (logFC)
# Both DI and DE upregulated (+logFC) in sulfidic populations.

DI_sig_DE_sig_upreg <- DI_sig_DE_sig %>% filter(logfc_peak>0 & logfc_rna>0)
print(paste0("Upregulated DI and DE: ", nrow(DI_sig_DE_sig_upreg)))
# write_csv(DI_sig_DE_sig_upreg, "analysis/00_jointAnalyses/new_DI_sig_DE_sig_upreg_01.15.24.csv")

# Summarize the annotations of both DI and DE upregulated (i.e., what region of the genome the peaks map to).

# Sum occurrences by genomic peak location
DI_sig_DE_sig_upreg_gen_loc <- as.data.frame(table(DI_sig_DE_sig_upreg$annotation))

# Change column header
names(DI_sig_DE_sig_upreg_gen_loc)[1] <- "annotation"
DI_sig_DE_sig_upreg_gen_loc

# Add additional columns
DI_sig_DE_sig_upreg_gen_loc$DE <- "1_DE_genes"
DI_sig_DE_sig_upreg_gen_loc$Direction <- "Upregulated"
DI_sig_DE_sig_upreg_gen_loc

# Add column of percentages
DI_sig_DE_sig_upreg_gen_loc_w_percents <- DI_sig_DE_sig_upreg_gen_loc %>% 
  mutate(Perc = `Freq` / sum(`Freq`)) %>% 
  arrange(Perc) %>% 
  mutate(Labels = scales::percent(Perc))
DI_sig_DE_sig_upreg_gen_loc_w_percents

# G-test (log likelihood ratio) of goodness of fit for both DI and DE upregulated. Do the annotations fit the distribution of all significantly DI peaks that have significant DE?
# Expected **percentages**
exp3_exon <- DI_sig_DE_sig_gen_loc_w_percents[DI_sig_DE_sig_gen_loc_w_percents$annotation=='Exon',3] %>% as.numeric() 
exp3_intergenic <- DI_sig_DE_sig_gen_loc_w_percents[DI_sig_gen_loc$annotation=='Intergenic',3] %>% as.numeric() 
exp3_intron <- DI_sig_DE_sig_gen_loc_w_percents[DI_sig_DE_sig_gen_loc_w_percents$annotation=='Intron',3] %>% as.numeric() 
exp3_prom <- DI_sig_DE_sig_gen_loc_w_percents[DI_sig_DE_sig_gen_loc_w_percents$annotation=='Promoter-TSS',3]%>% as.numeric() 
exp3_tts <- DI_sig_DE_sig_gen_loc_w_percents[DI_sig_DE_sig_gen_loc_w_percents$annotation=='TTS',3] %>% as.numeric() 

expected3_all_DI_peaks <- c(exp3_exon,exp3_intergenic,exp3_intron,exp3_prom,exp3_tts)
expected3_all_DI_peaks

# Observed **frequencies**
obs3_exon <- DI_sig_DE_sig_upreg_gen_loc_w_percents[DI_sig_DE_sig_upreg_gen_loc_w_percents$annotation=='Exon',2] %>% as.double() 
obs3_intergenic <- DI_sig_DE_sig_upreg_gen_loc_w_percents[DI_sig_DE_sig_upreg_gen_loc_w_percents$annotation=='Intergenic',2] %>% as.double() 
obs3_intron <- DI_sig_DE_sig_upreg_gen_loc_w_percents[DI_sig_DE_sig_upreg_gen_loc_w_percents$annotation=='Intron',2] %>% as.double() 
obs3_prom <- DI_sig_DE_sig_upreg_gen_loc_w_percents[DI_sig_DE_sig_upreg_gen_loc_w_percents$annotation=='Promoter-TSS',2] %>% as.double()
obs3_tts <- DI_sig_DE_sig_upreg_gen_loc_w_percents[DI_sig_DE_sig_upreg_gen_loc_w_percents$annotation=='TTS',2] %>% as.double() 

observed3_DI_sig_DE_sig <- c(obs3_exon,obs3_intergenic,obs3_intron,obs3_prom,obs3_tts)

GTest(x = observed3_DI_sig_DE_sig,
      p = expected3_all_DI_peaks,
      correct = "none")
# G = 1.6594, X-squared df = 4, p-value = 0.7981 #### NOT SIG


# Both DI and DE downregulated (-logFC) in sulfidic populations.

DI_sig_DE_sig_downreg <- DI_sig_DE_sig %>% filter(logfc_peak<0 & logfc_rna<0)
print(paste0("Downregulated DE and DI: ", nrow(DI_sig_DE_sig_downreg)))
# write.csv(x = DI_sig_DE_sig_downreg, file = "03_DI_sig_DE_sig_downreg.csv")

# Summarize the annotations of both DI and DE downregulated (i.e., what region of the genome the peaks map to).

# Sum occurrences by genomic peak location
DI_sig_DE_sig_downreg_gen_loc <- as.data.frame(table(DI_sig_DE_sig_downreg$annotation))
# Change column header
names(DI_sig_DE_sig_downreg_gen_loc)[1] <- "annotation"
DI_sig_DE_sig_downreg_gen_loc

# Add additional columns
DI_sig_DE_sig_downreg_gen_loc$DE <- "1_DE_genes"
DI_sig_DE_sig_downreg_gen_loc$Direction <- "Downregulated"
DI_sig_DE_sig_downreg_gen_loc

# Add column of percentages
DI_sig_DE_sig_downreg_gen_loc_w_percents <- DI_sig_DE_sig_downreg_gen_loc %>% 
  mutate(Perc = `Freq` / sum(`Freq`)) %>% 
  arrange(Perc) %>% 
  mutate(Labels = scales::percent(Perc))
DI_sig_DE_sig_downreg_gen_loc_w_percents

# G-test (log likelihood ratio) of goodness of fit for both DI and DE downregulated. Do the annotations fit the distribution of all significantly DI peaks that have significant DE?

# Expected **percentages**
expected4_all_DI_peaks <- expected3_all_DI_peaks

# Observed **frequencies**
obs4_exon <- DI_sig_DE_sig_downreg_gen_loc_w_percents[DI_sig_DE_sig_downreg_gen_loc_w_percents$annotation=='Exon',2] %>% as.double() 
obs4_intergenic <- DI_sig_DE_sig_downreg_gen_loc_w_percents[DI_sig_DE_sig_downreg_gen_loc_w_percents$annotation=='Intergenic',2] %>% as.double() 
obs4_intron <- DI_sig_DE_sig_downreg_gen_loc_w_percents[DI_sig_DE_sig_downreg_gen_loc_w_percents$annotation=='Intron',2] %>% as.double() 
obs4_prom <- DI_sig_DE_sig_downreg_gen_loc_w_percents[DI_sig_DE_sig_downreg_gen_loc_w_percents$annotation=='Promoter-TSS',2] %>% as.double()
obs4_tts <- DI_sig_DE_sig_downreg_gen_loc_w_percents[DI_sig_DE_sig_downreg_gen_loc_w_percents$annotation=='TTS',2] %>% as.double() 

observed4_DI_sig_DE_sig <- c(obs4_exon,obs4_intergenic,obs4_intron,obs4_prom,obs4_tts)

GTest(x = observed4_DI_sig_DE_sig,
      p = expected4_all_DI_peaks,
      correct = "none")
# G = 6.4474, X-squared df = 4, p-value = 0.1681 ## NOT SIG


# Divergent Peak and Gene expression.
# Combine divergent DI and DE in sulfidic populations.

DI_sig_DE_sig_divergent <- DI_sig_DE_sig %>% filter(logfc_peak * logfc_rna < 0)

# Summarize the annotations of DI and DE divergent (i.e., what region of the genome the peaks map to).

# Sum occurrences by genomic peak location
DI_sig_DE_sig_divergent_gen_loc <- as.data.frame(table(DI_sig_DE_sig_divergent$annotation))
# Change column header
names(DI_sig_DE_sig_divergent_gen_loc)[1] <- "annotation"
DI_sig_DE_sig_divergent_gen_loc

# Add additional columns
DI_sig_DE_sig_divergent_gen_loc$DE <- "1_DE_genes"
DI_sig_DE_sig_divergent_gen_loc$Direction <- "Divergent"
DI_sig_DE_sig_divergent_gen_loc

# Add column of percentages
DI_sig_DE_sig_divergent_gen_loc_w_percents <- DI_sig_DE_sig_divergent_gen_loc %>% 
  mutate(Perc = `Freq` / sum(`Freq`)) %>% 
  arrange(Perc) %>% 
  mutate(Labels = scales::percent(Perc))
DI_sig_DE_sig_divergent_gen_loc_w_percents

# G-test (log likelihood ratio) of goodness of fit for DI and DE divergent. Do the annotations fit the distribution of all significantly DI peaks that have significant DE?

expected5_all_DI_peaks <- expected3_all_DI_peaks

# Observed **frequencies**
obs5_exon <- DI_sig_DE_sig_divergent_gen_loc_w_percents[DI_sig_DE_sig_divergent_gen_loc_w_percents$annotation=='Exon',2] %>% as.double() 
obs5_intergenic <- DI_sig_DE_sig_divergent_gen_loc_w_percents[DI_sig_DE_sig_divergent_gen_loc_w_percents$annotation=='Intergenic',2] %>% as.double() 
obs5_intron <- DI_sig_DE_sig_divergent_gen_loc_w_percents[DI_sig_DE_sig_divergent_gen_loc_w_percents$annotation=='Intron',2] %>% as.double() 
obs5_prom <- DI_sig_DE_sig_divergent_gen_loc_w_percents[DI_sig_DE_sig_divergent_gen_loc_w_percents$annotation=='Promoter-TSS',2] %>% as.double()
obs5_tts <- DI_sig_DE_sig_divergent_gen_loc_w_percents[DI_sig_DE_sig_divergent_gen_loc_w_percents$annotation=='TTS',2] %>% as.double() 

observed5_DI_sig_DE_sig <- c(obs5_exon,obs5_intergenic,obs5_intron,obs5_prom,obs5_tts)

GTest(x = observed5_DI_sig_DE_sig,
      p = expected5_all_DI_peaks,
      correct = "none")
# G = 16.85, X-squared df = 4, p-value = 0.002067 ### SIG


# Generate a stacked bar chart of genomic locations (using a colorblind friendly palette).

# Combine data frames
gen_loc_w_percents2 <- rbind(DI_sig_DE_sig_upreg_gen_loc_w_percents, DI_sig_DE_sig_downreg_gen_loc_w_percents, DI_sig_DE_sig_divergent_gen_loc_w_percents)
gen_loc_w_percents2$DE <- factor(gen_loc_w_percents2$DE,levels=c('Divergent','Downregulated','Upregulated'))

# Color palette
colorBlind2   <- c("#E69F00", "#56B4E9", "#D55E00", "#009E73", "#CC79A7", "#F0E442")

# Plot
stacked_bar2 <- ggplot(data = gen_loc_w_percents2, aes(fill = annotation, x = reorder(Direction, desc(Direction)), y = Perc)) +
  # Stacked bar chart, use position = "stack" for frequency and position = "fill" for percentages
  geom_bar(position = "fill", stat = "identity", width = 0.8,show.legend = F) +
  #scale_fill_viridis(discrete = TRUE, option = "cividis") +
  scale_fill_manual(values = colorBlind) +
  scale_y_reverse() +
  xlab("") +
  ylab("Frequency") +
  theme(aspect.ratio = .1) +
  theme_classic(base_size = 14) +
  coord_flip()
stacked_bar2

library(patchwork)

stacked_bar + stacked_bar2

# Nodes and edges used to make an alluvial plot (below)
# Nodes and edges for significant DI but not significant DE.

# Add additional columns
DE <- rep("2_Not_DE", 5)
DI_sig_DE_nonsig_gen_loc2 <- cbind(DI_sig_DE_nonsig_gen_loc, DE)
DI_sig_DE_nonsig_gen_loc2$Direction <- "4_None"
DI_sig_DE_nonsig_gen_loc2


# Combine nodes and edges into a data frame that can be used to generate the alluvial plot.

# Bind nodes and edges
alluvial_df <- rbind(DI_sig_DE_sig_upreg_gen_loc, DI_sig_DE_sig_downreg_gen_loc, DI_sig_DE_sig_divergent_gen_loc, DI_sig_DE_nonsig_gen_loc2)

# Re-order columns
alluvial_df <- alluvial_df[, c(1, 3, 4, 2)]
alluvial_df

### Plot alluvial plot
alluvial_df_x3 <- rbind(alluvial_df, alluvial_df, alluvial_df)
alluvial_df_x3$id <- rep(c(1:20), times = 3)
alluvial_df_x3$x <- c(rep("Annotation", 20), rep("DE", 20), rep("Direction", 20))
alluvial_df_x3$y <- c(rep(c("Exon", "Intergenic", "Intron", "Promoter-TSS", "TTS"), 4),
                      rep("1_DE_genes", 15), rep("2_Not_DE", 5),
                      rep("1_Upregulated", 5), rep("2_Downregulated", 5), rep("3_Divergent", 5), rep("4_None", 5))              

ggplot(alluvial_df_x3, aes(x, id = id, split = y, value = Freq)) +
  geom_parallel_sets(aes(fill = annotation), alpha = 0.3, axis.width = 0.1,show.legend = F) +
  geom_parallel_sets_axes(aes(fill=y), axis.width = 0.1,show.legend = F) +
  geom_parallel_sets_labels(colour = 'black', size = 4,show.legend = F) + 
  theme_void()



