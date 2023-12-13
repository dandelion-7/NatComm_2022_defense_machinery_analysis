# This script follows script 9 to summarize the composition of metagenmoic samples after mapping reads with bowtie2 to assembled contigs annotated with mmseqs2
#---------------------------------------------------------------------------------------------------------------------
getwd()
library(dplyr)
library(ggplot2)
library(tidyr)
library(randomcoloR)
setwd('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts')

# loading taxonomy assignments of contigs
total_contigs <- read.table('../intermediates/15.mmseqs_taxonomy/taxonomy_new/mmseqs2_total_summary.txt', 
                      sep = '\t', header = T)

# pre-processing
mapping <- read.table('../intermediates/9.bowtie2/bam/total_stats.txt', 
                    header = F, sep = '\t')
sample_info <- as.data.frame(mapping[3])
sample_info <- sample_info %>% separate(V3, sep = '_', c('individual', 'time_point', 'sample'))
mapping <- cbind(mapping, sample_info)
colnames(mapping)[1:3] <- c('mapping_count', 'contig', 'sample_info')
mapping$contig <- paste(mapping$individual, mapping$contig, sep = '_')

contigs <- total_contigs[, 8:16]
contigs$contig <- paste(contigs$individual, contigs$contig_info, sep = '_')
contigs <- contigs[, -c(1, 2)]

total_stats <- left_join(mapping, contigs)
total_stats[, 7:13][is.na(total_stats[, 7:13])] <- 'unclassified'

total_stats$mapping_count <- as.numeric(total_stats$mapping_count)
total_stats$individual <- as.numeric(total_stats$individual)
total_stats$time_point <- as.numeric(total_stats$time_point)

sample_sum <- total_stats %>% group_by(sample_info) %>% 
  summarise(sample_count = sum(mapping_count))
total_stats <- left_join(total_stats, sample_sum)
write.table(total_stats, "~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/mmseqs2+bowtie2/mmseqs2+bowtie2_summarized_reports/20-5-1.bowtie2_summarized_total_report.txt", 
            sep = '\t', col.names = T, row.names = F, quote = F)
# save the total report.

# phylum level stats
phylum <- total_stats[, c(1, 2, 3, 4, 5, 8)]
phylum_sum <- phylum %>% group_by(individual, time_point, sample_info, phylum) %>% 
  summarise(phylum_count = sum(mapping_count))
phylum_sum <- left_join(phylum_sum, sample_sum)
phylum_sum <- phylum_sum %>% mutate(phylum_propotion = phylum_count/sample_count)
phylum_sum <- phylum_sum %>% mutate(phylum_type = case_when(
  phylum == 'unclassified' ~ 'unclassified', 
  phylum_propotion >= 0.01 ~ phylum,
  phylum_propotion < 0.01 ~ 'others'
))
n_distinct(phylum_sum$phylum_type)
phylum_sum$phylum_type <- factor(phylum_sum$phylum_type, 
                             levels = c('p__Actinomycetota', 'p__Bacillota', 
                                        'p__Bacteroidota', 'p__Euryarchaeota', 
                                        'p__Pseudomonadota', 'p__Thermodesulfobacteriota',
                                        'p__Uroviricota', 'p__Verrucomicrobiota', 
                                        'others', 'unclassified'))
#color_10 <- distinctColorPalette(8)
#color_10[9:10] <- c("#bdbdbd", "#737373")
color_10 <- c("#8CAED5", "#D97B5F", "#D7D762", "#B758D6", "#7EE075", 
              "#D8D5BC", "#D995CA", "#84E0D0", "#bdbdbd", "#737373")
phylum_sum %>% ggplot(aes(x = time_point, y = phylum_propotion)) + 
  geom_col(aes(fill = phylum_type)) + 
  facet_wrap(. ~ individual, nrow = 2) + 
  scale_fill_manual(values = color_10) + 
  scale_x_continuous(limits = c(0.5, 10.5), 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/mmseqs2+bowtie2/mmseqs2+bowtie2_figures/20-5-1.phylum_level_composition
write.table(phylum_sum, '~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/mmseqs2+bowtie2/mmseqs2+bowtie2_summarized_reports/20-5-2.bowtie2_phylum_report.txt', 
            sep = '\t', col.names = T, row.names = F, quote = F)
# save the phylum-level report.

# genus level stats
genus <- total_stats[, c(1, 2, 3, 4, 5, 12)]
genus_sum <- genus %>% group_by(individual, time_point, sample_info, genus) %>% 
  summarise(genus_count = sum(mapping_count))
genus_sum <- left_join(genus_sum, sample_sum)
genus_sum <- genus_sum %>% mutate(genus_proportion = genus_count/sample_count)
genus_sum <- genus_sum %>% mutate(genus_type = case_when(
  genus == 'unclassified' ~ 'unclassified', 
  genus_proportion >= 0.01 ~ genus,
  genus_proportion < 0.01 ~ 'others'
))
n_distinct(genus_sum$genus_type)
#color_43 <- distinctColorPalette(41)
#color_43[42:43] <- c("#bdbdbd", "#737373")
color_43 <- c("#F135CC", "#E1E4CF", "#A7AFDF", "#E596E2", "#86EB2D", "#DBE5A5", 
              "#DFA4A0", "#95D9E6", "#A681DE", "#E7509A", "#5563DA", "#B9EAC9", 
              "#67E6C4", "#C976E8", "#5DE55F", "#E666DC", "#D7D569", "#D36C51", 
              "#8E9DEC", "#D5E3F0", "#63E79C", "#9A3C8F", "#9DE674", "#68A594",
              "#565DA2", "#DAC499", "#ACE599", "#E43E53", "#5EA9DA", "#DCC5D9",
              "#CD2FED", "#894BDB", "#C89A61", "#69E7E7", "#E7AFDA", "#DDEA4E", 
              "#6B9B4F", "#837F84", "#633BEE", "#D47497", "#E6A637", "#bdbdbd", 
              "#737373")
genus_sum %>% ggplot(aes(x = time_point, y = genus_proportion)) + 
  geom_col(aes(fill = genus_type)) + 
  facet_wrap(. ~ individual, nrow = 2) + 
  scale_fill_manual(values = color_43) + 
  scale_x_continuous(limits = c(0.5, 10.5), 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/mmseqs2+bowtie2/mmseqs2+bowtie2_figures/20-5-2.genus_level_composition

genus_sum %>% filter(genus_type != 'unclassified') %>% 
  ggplot(aes(x = time_point, y = genus_proportion)) + 
  geom_col(aes(fill = genus_type)) + 
  facet_wrap(. ~ individual, nrow = 2) + 
  scale_fill_manual(values = color_43[1:42]) + 
  scale_x_continuous(limits = c(0.5, 10.5), 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/mmseqs2+bowtie2/mmseqs2+bowtie2_figures/20-5-3.genus_level_composition_no_unclassified
write.table(genus_sum, '~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/mmseqs2+bowtie2/mmseqs2+bowtie2_summarized_reports/20-5-3.bowtie2_genus_report.txt', 
            sep = '\t', col.names = T, row.names = F, quote = F)
# save the genus-level report.
