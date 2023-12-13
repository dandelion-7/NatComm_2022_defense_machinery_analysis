# This script is for summarizing the kmers of stable repeats from each individual, calculating the BC distance, and perform NMDS to visualize the differences of stable CRISPR repeats among individuals.
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

getwd()
setwd('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts')

library(vegan)
library(ggplot2)
library(tidyverse)

kmer_list <- read.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/11.clustered_CRISPR_stats/11_2.clustered_repeats_kmer_counting/stable_repeats_5_mer_summary.txt',
                        sep = '\t', header = FALSE)

colnames(kmer_list) <- c('kmer', 'count', 'individual', 'timepoint')
kmer_list <- kmer_list %>% mutate(sample = paste(individual, timepoint, sep = '_'))

kmer_env <- kmer_list %>% group_by(individual, timepoint, sample) %>% 
  summarise(count = 1)

kmer_list <- kmer_list[c(-3, -4)]

kmer_summary <- kmer_list %>% group_by(kmer, sample) %>% 
  summarise(count = n())

kmer_matrix <- spread(kmer_list, key = "sample", value = "count")
kmer_matrix <- t(kmer_matrix)
colnames(kmer_matrix) <- kmer_matrix[1, ]
kmer_matrix <- kmer_matrix[-1, ]
kmer_matrix <- as.data.frame(kmer_matrix)
kmer_matrix[is.na(kmer_matrix)] <- 0

type(kmer_matrix)
kmer_matrix <- apply(kmer_matrix, 2, as.numeric)

distance <- vegdist(kmer_matrix, method = 'bray')
nmds <- metaMDS(distance, k = 2)
kmer_stress <- nmds$stress
kmer_df <- as.data.frame(nmds$points)
kmer_df <- cbind(kmer_df, kmer_env)

fill_colors <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')
kmer_df$individual <- as.character(kmer_df$individual)
kmer_df$individual <- factor(kmer_df$individual, levels = c('1', '2', '3', '4', 
                                                            '5', '6', '7', '8',
                                                            '9', '10'))
kmer_df <- kmer_df[order(kmer_df$individual),]

kmer_df %>% ggplot(aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(fill = individual), 
             shape = 21, 
             color = 'black', 
             size = 3, 
             alpha = 0.8) + 
  stat_ellipse(aes(color = individual), level = 0.95, size = 0.8) + 
  scale_fill_manual(values = fill_colors[c(1:10)]) + 
  scale_color_manual(values = fill_colors[c(1:10)]) + 
  guides(color = 'none') + 
  labs(x = 'NMDS1', y = 'NMDS2', fill = 'Individual')
# clustered_repeats_kmer_profiling
