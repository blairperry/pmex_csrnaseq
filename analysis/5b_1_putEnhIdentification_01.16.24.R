

library(tidyverse)
library(cols4all)

sulf <- c('MX31','MX32','MX46','MX48','MX76','MX77')

peak_ints <- read_tsv('data/csRNA/07_homer_mergepeaks/allintersects.pmex.txt',
                      col_names = c('chr_mergedpeak','start_mergedpeak','end_mergedpeak','id_mergedpeak',
                                    'score_mergedpeak','strand_mergedpeak','file_peak','chr_peak','start_peak',
                                    'end_peak','id_peak','score_peak','strand_peak','annotation_peak',
                                    'nearestGene_peak','nearestGeneDist_peak','proximity_peak','stability_peak'))

peak_simple <- peak_ints %>%
  select(id_mergedpeak,file_peak) %>% 
  unique() %>% 
  group_by(id_mergedpeak) %>% 
  add_count()

peak_hist <-  peak_simple %>% 
  group_by(n) %>% 
  add_count() %>% 
  select(n,nn) %>% 
  unique()

ggplot(peak_hist,aes(x=n,y=nn)) +
  geom_bar(stat='identity') +
  labs(x='# Samples in Intersection',y='# Merged Peaks') +
  scale_x_continuous(breaks = c(1:11)) +
  theme_linedraw(base_size = 16) 

sample_hist <- peak_simple %>% 
  ungroup() %>% 
  select(1,2) %>% 
  group_by(file_peak) %>% 
  mutate(file_peak = str_split_fixed(file_peak,'[/]',4)[,4] %>% str_remove_all('.tss.txt.bed')) %>% 
  add_count() %>% 
  select(file_peak,n) %>% 
  unique() %>% 
  mutate(condition = ifelse(file_peak %in% sulf,'sulfidic','non-sulfidic'))

ggplot(sample_hist,aes(x=file_peak,y=n,fill=condition)) +
  geom_bar(stat='identity') +
  labs(x='Sample',y='# Peaks') +
  theme_linedraw() 

peak_simple %>% 
  ungroup() %>% 
  select(2,3) %>% 
  mutate(file_peak = str_split_fixed(file_peak,'[/]',4)[,4] %>% str_remove_all('.tss.txt.bed')) %>% 
  group_by(file_peak,n) %>% 
  tally() %>% 
  ggplot(aes(x=n,y=nn,fill=file_peak)) +
  geom_bar(stat='identity',position = 'fill') +
  scale_x_continuous(breaks = c(1:12)) +
  # scale_fill_discrete_c4a_cat("carto.safe") +
  theme_linedraw()

peak_simple %>% 
  ungroup() %>% 
  select(2,3) %>% 
  mutate(file_peak = str_split_fixed(file_peak,'[/]',4)[,4] %>% str_remove_all('.tss.txt.bed')) %>% 
  group_by(file_peak,n) %>% 
  tally() %>% 
  ggplot(aes(x=n,y=nn,fill=file_peak)) +
  geom_bar(stat='identity') +
  scale_x_continuous(breaks = c(1:12)) +
  # scale_fill_discrete_c4a_cat("carto.safe") +
  theme_linedraw()


peak_simple_out <- peak_simple %>% select(mergedpeak_id=1,intersect_id=2,num_intersect_samples=3) %>% 
  mutate(intersect_id = str_split_fixed(intersect_id,'[/]',4)[,4] %>% str_remove_all('.tss.txt.bed'))
# write_tsv(peak_simple_out,'analysis/1_peakExploration/peakIntersectCounts_01.12.24.txt')

## Making table of peak annotations

peak_ints.types <- peak_ints %>% 
  select(id_mergedpeak,annotation_peak,nearestGene_peak,proximity_peak,stability_peak) %>% 
  unique() %>% 
  group_by(id_mergedpeak) %>% add_count(name = 'annotation_count')

peak_ints.types %>% filter(annotation_count==1) # 54,146 merged peaks have consistent annotations across all involved sub-peaks 
(peak_ints.types %>% filter(annotation_count==1) %>% nrow()) / (peak_ints %>% select(id_mergedpeak) %>% unique() %>% nrow()) # equates to ~81%

# write_tsv(peak_ints.types,'analysis/1_peakExploration/peakIntersectAnnotations_01.12.24.bed')


# Make table of putative enhancer peaks -----------------------------------

# first filter by peaks with strand-aware annotations of distal, unstable, and other
put_enh <- peak_ints %>% 
  filter(proximity_peak=='distal' & stability_peak=='unstable' & annotation_peak=='other') %>% 
  filter(abs(nearestGeneDist_peak) >= 500) %>% # peaks must be at least 500 bp away from nearest gene (i.e., outside likely promoter)
  select(chr_mergedpeak,start_mergedpeak,end_mergedpeak,id_mergedpeak) %>% 
  unique() 
 
nrow(put_enh) # 33,799 putative enhancer peaks

# write_tsv(put_enh,'./analysis/1_PeakExploration/putativeEnhancerPeaks_pmex_01.12.24.bed',col_names = F)

