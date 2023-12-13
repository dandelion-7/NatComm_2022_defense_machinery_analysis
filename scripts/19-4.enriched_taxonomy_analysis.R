# This script is for comparing the differential contributions of taxonomies to the total microbiota and the CRISPR repertoire, following script 19-2. But the paired t-test of the few varied taxonomies are not done yet.
#-------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(randomcoloR)

getwd()
setwd('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/')

total <- read.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/19.bowtie2_reads_contig/figures/19-3-1.crispr_reads_taxonomy_stat.txt', 
                    sep = '\t', header = T)
CRISPR_genus_stat <- total %>% group_by(individual, time_point, sample, 
                                        sample_info, genus) %>% 
  summarise(genus_count = n())
sample_read_sum <- total %>% group_by(sample_info) %>% 
  summarise(read_count = n())
CRISPR_genus_stat <- left_join(CRISPR_genus_stat, sample_read_sum, by = 'sample_info')
CRISPR_genus_stat <- CRISPR_genus_stat %>% mutate(genus_proportion = genus_count/read_count)
CRISPR_genus_stat$object <- "CRISPR_reads"

bowtie2_genus <- read.table("~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/mmseqs2+bowtie2/mmseqs2+bowtie2_summarized_reports/20-5-3.bowtie2_genus_report.txt", 
                            header = T, sep = '\t')
bowtie2_genus$object <- "microbiome_reads"
bowtie2_genus <- bowtie2_genus[, c(1, 2, 3, 4, 7, 8, 9)]

CRISPR_genus_average <- CRISPR_genus_stat %>% group_by(genus) %>% 
  summarize(average_CRISPR_proportion = sum(genus_proportion)/100)

bowtie2_genus_average <- bowtie2_genus %>% group_by(genus) %>% 
  summarize(average_bowtie2_proportion = sum(genus_proportion)/100)

average_genus_comparing <- left_join(CRISPR_genus_average, bowtie2_genus_average) %>% drop_na()
average_genus_comparing <- average_genus_comparing %>% mutate(CRISPR_enrichment = average_CRISPR_proportion/average_bowtie2_proportion)
average_genus_comparing %>% filter(genus != 'unclassified') %>%
  ggplot(aes(x = average_bowtie2_proportion, y = average_CRISPR_proportion)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  geom_vline(xintercept = 0.005, linetype = 'dashed') + 
  geom_hline(yintercept = 0.005, linetype = 'dashed') +
  scale_x_continuous(limits = c(-0.00001, 0.09), breaks = c(0, 0.005, 0.03, 0.06, 0.09)) + 
  scale_y_continuous(limits = c(-0.00001, 0.09), breaks = c(0, 0.005, 0.03, 0.06, 0.09)) + 
  scale_fill_manual(values = color_20)
#.../19.bowtie2_reads_contig/figures/19-4-1.genus_average_proportion_comparing

#color_20 <- distinctColorPalette(20)
color_20 <- c("#D0E5DA", "#D3B14D", "#81EABC", "#726BD9", "#DB4390",
              "#E38FBA", "#846785", "#71B5D5", "#74E2E3", "#D583E1", 
              "#6E997C", "#DD6E51", "#71EA4A", "#92A2E1", "#DAC1DF", 
              "#CBE556", "#CCE199", "#DBB196", "#BF43E0", "#71DC7C")
average_genus_comparing %>% filter(genus != 'unclassified') %>% 
  filter(average_bowtie2_proportion >= 0.005 | average_CRISPR_proportion >= 0.005) %>% 
  ggplot(aes(x = reorder(genus, -CRISPR_enrichment), y = CRISPR_enrichment)) +
  geom_col(aes(fill = genus)) + 
  geom_hline(yintercept = 1, linetype = 'dashed') + 
  scale_fill_manual(values = color_20) + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, size = 8)) + 
  theme(legend.position = 'None')
#.../19.bowtie2_reads_contig/figures/19-4-2.proportion_variation_of_main_taxonomies

combined_genus <- rbind(bowtie2_genus[, c(3, 4, 5, 7)], CRISPR_genus_stat[, c(4, 5, 8, 9)])
combined_genus <- left_join(combined_genus, average_genus_comparing, by = 'genus')
combined_genus$object <- factor(combined_genus$object, levels = c('microbiome_reads', 'CRISPR_reads'))
combined_genus %>% filter(genus != 'unclassified') %>% 
  filter(average_CRISPR_proportion >= 0.005 | average_bowtie2_proportion >= 0.005) %>% 
  ggplot(aes(x = object, y = genus_proportion)) + 
  geom_boxplot() + 
  geom_point(shape = 21, fill = 'grey', size = 2, position = position_dodge(0.3), alpha = 0.8) + 
  geom_line(aes(group = sample_info, color = genus), alpha = 0.8) + 
  facet_wrap(.~ genus, nrow = 4) + 
  scale_color_manual(values = color_20)+ 
  theme(legend.position = 'None') + 
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1, size = 8))
#.../19.bowtie2_reads_contig/figures/19-4-3.proportion_comparison_of_main_taxonomies

# perform paired-T test on the few differed taxonomies: Bacteroides/Escherichia/Ligilactobacillus/Megamonas/Ruminococcus
microbiome <- combined_genus[combined_genus$object == 'microbiome_reads' & combined_genus$genus == 'g__Ligilactobacillus',] %>% drop_na()
crispr <- combined_genus[combined_genus$object == 'CRISPR_reads' & combined_genus$genus == 'g__Ligilactobacillus',] %>% drop_na()
microbiome_crispr <- full_join(crispr[, c(1,3)], microbiome[, c(1,3)], by = 'sample_info')
microbiome_crispr[is.na(microbiome_crispr)] <- 0
t.test(microbiome$genus_proportion, crispr$genus_proportion, paired = F)
