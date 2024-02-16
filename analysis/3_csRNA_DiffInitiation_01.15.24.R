
library(tidyverse)
library(edgeR)
library(colorspace)

sample_info <- readxl::read_xlsx('analysis/sampleInfo_simple_01.15.24.xlsx',col_names = c('sample_id','species','drainage','ecotype'))

all_peak_data <- read_tsv(file = 'data/csRNA/08_homer_annotatepeaks/rawcounts.txt') %>% 
  select(peak_id=1,everything())

raw_data <- all_peak_data %>% 
  select(peak_id=1,contains('normalization')) %>% 
  pivot_longer(-1,names_to = 'sample',values_to = 'count') %>% 
  mutate(sample = str_split_fixed(sample,'[_]',8)[,7]) %>% 
  filter(sample %in% sample_info$sample_id) %>% 
  left_join(sample_info,by=c('sample'='sample_id')) %>% 
  mutate(sample_info = paste(sample,species,drainage,ecotype,sep = ':')) %>% 
  select(peak_id,count,sample_info) %>% 
  pivot_wider(names_from = sample_info,values_from = count)

group = str_split_fixed(names(raw_data[,-1]),'[:]',4)[,4]
drainage = str_split_fixed(names(raw_data[,-1]),'[:]',4)[,3]
species = str_split_fixed(names(raw_data[,-1]),'[:]',4)[,2]

raw_counts <- as.data.frame(raw_data) %>% column_to_rownames('peak_id')

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
# Samples must have a cpm of at least 2.104157 (~5 counts in smallest library) and be present in at least 5 samples
keep <- rowSums(cpm(y) >= 2.104157) >= 5
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
res.df <- data.frame(res) %>% rownames_to_column('peak_id')

res.df %>% filter(FDR<0.05) %>% nrow()
res.df.sig <- res.df %>% filter(FDR<0.05) 

res.df_fullInfo <- res.df %>% left_join(all_peak_data,by='peak_id')
res.df.sig_fullInfo <- res.df.sig %>% left_join(all_peak_data,by='peak_id')
nrow(res.df.sig_fullInfo)


# write_tsv(res.df_fullInfo,'analysis/0_DEanalysis_new/csRNA_DEResults_All_01.15.24.tsv')
# write_tsv(res.df.sig_fullInfo,'analysis/0_DEanalysis_new/csRNA_DEResults_Sig_01.15.24.tsv')

# Add human annotation info and format for supp table --------------

human_annot <- readxl::read_xlsx('/mec14360-sup-0002-tables2.xlsx') %>% 
  janitor::clean_names() %>% 
  select(gene_name, human_id = 3, protein_annotations)

res.df.sig_fullInfo_humanInfo <- res.df.sig_fullInfo %>% left_join(human_annot,by=c('Gene Name'='gene_name')) %>% 
  mutate(human_homolog = str_split_fixed(human_id,'[|]',3)[,3] %>% str_remove_all('_HUMAN'))

names(res.df.sig_fullInfo_humanInfo)

supp_table_s5 <- res.df.sig_fullInfo_humanInfo %>% 
  arrange(-logFC) %>% 
  mutate(human_homolog = ifelse(human_homolog=='',NA,human_homolog)) %>% 
  select(`csRNA Peak ID` = peak_id,
         Chromosome = Chr,
         Start,End,Strand,
         Annotation,
         `Nearest Gene` = `Gene Name`,
         `Distance to TSS (bp)` = `Distance to TSS`,
         `Log2 Fold Change`=logFC,
         `p-value`=PValue,
         FDR,
         `Human Homologue` = human_homolog,
         `Human Protein Annotation`=protein_annotations
         )
# write_tsv(supp_table_s5,'analysis/0_DEanalysis_new/supp_table_s5_01.30.24.tsv')



###

cpm_csrna <- cpm(filtered_y, normalized.lib.sizes=TRUE, log=F)
cpm_csrna <- as.data.frame(cpm_csrna[,])

logcpm_csrna <- cpm(filtered_y, normalized.lib.sizes=TRUE, log=T,prior.count = 1)
logcpm_csrna <- as.data.frame(logcpm_csrna[,])

cpm_counts.sig <- logcpm_csrna %>% rownames_to_column('peak_id') %>% 
  filter(peak_id %in% res.df.sig$peak_id) %>% 
  column_to_rownames('peak_id')

cpm_counts.sig_top5k <- cpm_csrna %>% 
  mutate(sum=rowSums(.[,])) %>% 
  top_n(5000,wt = sum) %>% 
  select(-sum)

# Top 5k
pheatmap::pheatmap(log10(cpm_counts.sig_top5k+1),scale='row',
                   # cutree_cols = 2,
                   # cutree_rows = 2,
                   treeheight_row = 0,
                   show_rownames = F,
                   color = diverging_hcl(n = 5,palette = 'Blue-Red 2')
                   )

# DE
pheatmap::pheatmap(cpm_counts.sig,scale='row',
                   cutree_cols =2,
                   treeheight_row = 0,
                   show_rownames = F,
                   # gaps_col = c(2,4,6),
                   cutree_rows = 2,
                   color = diverging_hcl(n = 5,palette = 'Blue-Red 2'))


# Test overlap with candidate gene set ------------------------------------

# Merge significant peaks with H2S-related candidate list. The Sulfide Detox/Response Gene Set is relevant for this study. 
# The NuclearRef and Broughton Gene Sets are reference sets.

# Read in candidate list
cand_list <- read.csv('OXPHOS_Reference_and_Detox_Gene_IDS.csv', header = 1)

# Subset to only include Sulfide Detox/Response Gene Set
subset_cand <- subset(cand_list, Gene.Set == "Sulfide Detox/Response")

# Change name of sulfide:quinone oxidoreductase from "sqor" to "sqrdl" in Gene.name.from.accession column
subset_cand[subset_cand == "sqor"] <- "sqrdl"

# Merge H2S candidate list with DI peaks
res.df.sig_fullInfo_cand <- res.df.sig_fullInfo %>% 
  janitor::clean_names() %>% 
  filter(gene_name %in% subset_cand$Gene.name.from.accession)

# write_tsv(res.df.sig_fullInfo_cand,'analysis/0_DEanalysis_new/csRNA_DEresults_CandGenes_01.26.24.tsv')

focal <- res.df.sig_fullInfo_cand %>% filter(gene_name %in% c('sqrdl','LOC106917690','LOC106917689'))

# Determine if the number of differentially expressed peaks in the Sulfide Detox/Response Gene Set is significant using a Fisher's Exact Test.

#2 x 2 contingency table (Set = Sulfide Detox/Response Gene Set):

### DI peaks

di_cand <- nrow(res.df.sig_fullInfo_cand)
di_notcand <- nrow(res.df.sig_fullInfo) - nrow(res.df.sig_fullInfo_cand)

notDI_cand <- nrow(subset_cand) - nrow(res.df.sig_fullInfo_cand)
notDI_notcand <- res.df_fullInfo %>% 
  janitor::clean_names() %>% 
  filter(!(gene_name %in% subset_cand$Gene.name.from.accession)) %>% 
  filter(!(gene_name %in% res.df.sig_fullInfo$gene_name)) %>% 
  nrow()

# Make the 2 x 2 contingency table as a matrix
table_up <- matrix(c(di_cand, di_notcand, notDI_cand, notDI_notcand),
                   nrow = 2,
                   dimnames = list(c('in_set', 'not_in_set'),
                                   c('target_up', 'not_target_up')))
# Fisher's Exact Test for count data
fisher.test(table_up, alternative = "two.sided") # Sig, p=0.03229



####