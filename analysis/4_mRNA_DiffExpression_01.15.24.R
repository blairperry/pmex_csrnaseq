
library(tidyverse)
library(edgeR)

sample_info <- readxl::read_xlsx('analysis/sampleInfo_simple_01.15.24.xlsx',col_names = c('sample_id','species','drainage','ecotype'))

wild_data <- read_csv(file = 'data_fromKerryAnalysis/mRNA/raw_counts/gene_count_matrix-2019-10-01_no_STRG.csv') %>% 
  pivot_longer(-1,names_to = 'sample',values_to = 'count') %>% 
  filter(sample %in% sample_info$sample_id) %>% 
  left_join(sample_info,by=c('sample'='sample_id')) %>% 
  mutate(sample_info = paste(sample,species,drainage,ecotype,sep = ':')) %>% 
  select(gene_id,count,sample_info) %>% 
  pivot_wider(names_from = sample_info,values_from = count)

group = str_split_fixed(names(wild_data[,-1]),'[:]',4)[,4]
drainage = str_split_fixed(names(wild_data[,-1]),'[:]',4)[,3]
species = str_split_fixed(names(wild_data[,-1]),'[:]',4)[,2]

raw_counts <- as.data.frame(wild_data) %>% column_to_rownames('gene_id')

y <- DGEList(counts=raw_counts, group=group)

# Sort samples by library size
y_sorted <- y$samples[order(y$samples$lib.size) , ]
y_sorted
# Select lowest library size
lowest_library_size <- y_sorted$lib.size[1]
# Print lowest library size
cat("Smallest wild library size:", lowest_library_size)
# Calculate minimum counts per million (cpm) to have ~5 counts in the smallest library, as per edgeR manual recommendations
# (5)/(lowest_library_size)=(x_wild)/(1,000,000)
# Solve for x_wild
x <- (5*1000000)/lowest_library_size
# Print minimum cpm
cat("\n")
cat("Minimum wild cpm:", x)
# Samples must have a cpm of at least 1.144 (~5 counts in smallest library) and be present in at least 5 samples
keep <- rowSums(cpm(y) >= 1.144) >= 5
filtered_y <- y[keep, , keep.lib.sizes=FALSE]
# Library sizes post-filtering
filtered_y$samples

# Normalize for RNA composition
filtered_y <- calcNormFactors(filtered_y)


# Create design matrix.
design <- model.matrix(~group + drainage + species)

filtered_y <- estimateDisp(filtered_y, design)
BCV <- plotBCV(filtered_y)
fit <- glmQLFit(filtered_y, design)

lrt <- glmLRT(fit, coef = 2)
res <- topTags(lrt,n = Inf)
res.df <- data.frame(res) %>% rownames_to_column('gene_id')

res.df %>% filter(FDR<0.05) %>% nrow()
res.df.sig <- res.df %>% filter(FDR<0.05) 


# Add human annotation ----------------------------------------------------

pmex_name_convert <- read_csv('data_fromKerryAnalysis/reference/PmexGeneNameMatching.csv',skip = 1)
human_annot <- readxl::read_xlsx('data_fromKerryAnalysis/reference/mec14360-sup-0002-tables2.xlsx') %>% janitor::clean_names() %>% 
  select(gene_id,gene_name,subject_sequence_id,protein_annotations)

res.df_fullInfo <- res.df %>% 
  left_join(pmex_name_convert,by=c('gene_id'='gene.ID')) %>% 
  select(-B2Gname,-rnaID) %>% 
  left_join(human_annot,by='gene_id') %>% 
  select(gene_id,gene.name,human_id=subject_sequence_id,protein_annotations,everything(),-gene_name) %>% 
  dplyr::rename(gene_name=gene.name) %>% 
  unique()

res.df_fullInfo %>% group_by(gene_id) %>% add_count() %>% arrange(-n)

res.df.sig_fullInfo <- res.df.sig %>% 
  left_join(pmex_name_convert,by=c('gene_id'='gene.ID')) %>% 
  select(-B2Gname,-rnaID) %>% 
  left_join(human_annot,by='gene_id') %>% 
  select(gene_id,gene.name,human_id=subject_sequence_id,protein_annotations,everything(),-gene_name) %>% 
  unique()


write_tsv(res.df_fullInfo,'analysis/0_DEanalysis_new/mRNA_DEResults_All_01.15.24.tsv')
write_tsv(res.df.sig_fullInfo,'analysis/0_DEanalysis_new/mRNA_DEResults_Sig_01.15.24.tsv')



# Reformat results for Supp Table S8 --------------------------------------
names(res.df.sig_fullInfo)
supp_table_s8 <- res.df.sig_fullInfo %>% 
  mutate(human_homologue = str_split_fixed(human_id,'[|]',3)[,3]) %>% 
  mutate(human_homologue = ifelse(is.na(protein_annotations),NA,human_homologue)) %>% 
  select(`Gene Name`=gene.name,
         `Human Homologue`= human_homologue,
         `Human Protein Annotation`=protein_annotations,
         `Log2 Fold Change`=logFC,
         `p-value`=PValue,
         FDR)
write_tsv(supp_table_s8,'analysis/0_DEanalysis_new/supp_table_s8_01.30.24.tsv')


#

boxplot(as.numeric(raw_counts["gene11852", ]) ~ group)


# Compare to previous results ---------------------------------------------

de.old <- read_csv('data_fromKerryAnalysis/mRNA/edgeR_wgcna/5_habitat_qlf_S_vs_NS_WITH_ANNOTATIONS.csv') %>% 
  janitor::clean_names() %>% 
  select(-x1) %>% 
  filter(fdr<0.05)

length(row.names(res.df.sig) %in% de.old$gene_id) / length(row.names(res.df.sig))
length(de.old$gene_id)


####