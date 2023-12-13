# This script follows script 20-3/4/5, to combine the taxonomic results of Kraken2/Metaphlan4/MMseqs2+bowtie2, to compare if the MMseqs2+bowtie2 is comparable with the two classical softwares.
#-----------------------------------------------------------------------------------------------------------------
getwd()
setwd('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts')
library(ggplot2)
library(dplyr)
library(tidyr)
library(randomcoloR)

# kraken_phylum
kraken_phylum <- read.table("~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/kraken/summarized_reports/20-3-2.kraken_phylum_reports.txt",
                            header = T, sep = '\t')
kraken_phylum$method <- 'Kraken2'
kraken_phylum <- kraken_phylum[, c(7, 8, 9, 15, 16, 17, 18)]
# kraken_genus
kraken_genus <- read.table("~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/kraken/summarized_reports/20-3-3.kraken_genus_reports.txt", 
                           header = T, sep = '\t')
kraken_genus$method <- 'Kraken2'
kraken_genus <- kraken_genus[, c(7, 8, 9, 15, 16, 17, 18)]

# metaphlan_phylum
mpa_phylum <- read.table("~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/metaphlan/summarized_reports/20-4-2.mpa_phylum_report.txt", 
                         header = T, sep = '\t')
mpa_phylum$method <- "Metaphlan4"

mpa_phylum <- mpa_phylum[, c(4, 12, 13, 22, 3, 23, 24)]
# metaphlan_genus
mpa_genus <- read.table("~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/metaphlan/summarized_reports/20-4-3.mpa_genus_report.txt", 
                         header = T, sep = '\t')
mpa_genus$method <- "Metaphlan4"
mpa_genus <- mpa_genus[, c(4, 12, 13, 22, 3, 23, 24)]

# mmseqs2+bowtie2_phylum
bowtie2_phylum <- read.table("~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/mmseqs2+bowtie2/mmseqs2+bowtie2_summarized_reports/20-5-2.bowtie2_phylum_report.txt", 
                             header = T, sep = '\t')
bowtie2_phylum$method <- "mmseqs2+bowtie2"
bowtie2_phylum <- bowtie2_phylum[, c(3, 1, 2, 4, 7, 8, 9)]
# mmseqs2+bowtie2_genus
bowtie2_genus <- read.table("~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/mmseqs2+bowtie2/mmseqs2+bowtie2_summarized_reports/20-5-3.bowtie2_genus_report.txt", 
                             header = T, sep = '\t')
bowtie2_genus$method <- "mmseqs2+bowtie2"
bowtie2_genus <- bowtie2_genus[, c(3, 1, 2, 4, 7, 8, 9)]

# unify the column names
colnames(mpa_phylum) <- colnames(kraken_phylum)
colnames(bowtie2_phylum) <- colnames(kraken_phylum)
colnames(mpa_genus) <- colnames(kraken_genus)
colnames(bowtie2_genus) <- colnames(kraken_genus)

# convert the proportion of metaphlan from percentage to proportion.
mpa_phylum <- mpa_phylum %>% mutate(phylum_proportion = phylum_proportion/100)
mpa_genus <- mpa_genus %>% mutate(genus_proportion = genus_proportion/100)

# phylum-level comparing
phylum <- rbind(kraken_phylum, mpa_phylum, bowtie2_phylum)
n_distinct(phylum$phylum_type)
phylum$phylum_type <- factor(phylum$phylum_type, levels = c(
  "p__Actinomycetota", "p__Bacillota", "p__Bacteroidota", "p__Euryarchaeota",
  "p__Pseudomonadota", "p__Thermodesulfobacteriota", "p__Uroviricota", 
  "p__Verrucomicrobiota", "others", "unclassified"
))
color_10 <- c("#8CAED5", "#D97B5F", "#D7D762", "#B758D6", "#7EE075", 
              "#D8D5BC", "#D995CA", "#84E0D0", "#bdbdbd", "#737373")
phylum %>% ggplot(aes(x = time_point, y = phylum_proportion)) + 
  geom_col(aes(fill = phylum_type)) + 
  facet_grid(method ~ individual) +
  scale_fill_manual(values = color_10) + 
  scale_x_continuous(limits = c(0.5, 10.5), 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) #20-6-1.combined_phylum_level_comparison

# genus-level comparing
genus <- rbind(kraken_genus, mpa_genus, bowtie2_genus)
n_distinct(genus$genus_type)
#color_59 <- distinctColorPalette(57)
#color_59[58:59] <- c("#bdbdbd", "#737373")
genus %>% ggplot(aes(x = time_point, y = genus_proportion)) + 
  geom_col(aes(fill = genus_type)) + 
  facet_grid(method ~ individual) + 
  scale_fill_manual(values = color_59) + 
  scale_x_continuous(limits = c(0.5, 10.5), 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
  guides(fill = guide_legend(ncol = 1, )) + 
  theme(legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size = 10)) #20-6-2.combined_genus_level_comparison
