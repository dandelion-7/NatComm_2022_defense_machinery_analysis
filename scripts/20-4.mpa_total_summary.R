# This script is for summarizing the taxonomic annotation results from Metaphlan4, following script 20-2.
#--------------------------------------------------------------------------------------------
getwd()
library(dplyr)
library(tidyr)
library(ggplot2)
library(randomcoloR)
setwd('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts')

# pre-processing
total <- read.table('../intermediates/20.taxonomy_profiling/metaphlan/mpa/total_output.txt', 
                    sep = '\t', header = F)
taxonomy <- as.data.frame(total$V1)
taxonomy <- taxonomy %>% separate(`total$V1`, sep = ",", 
                                  c('kingdom', 'phylum', 'class', 'order', 
                                    'family', 'genus', 'species'))
taxid <- as.data.frame(total$V2)
taxid <- taxid %>% separate(`total$V2`, sep = ',',
                            c('kingdom_id', 'phylum_id', 'class_id', 'order_id', 
                              'family_id', 'genus_id', 'species_id'))
taxid[is.na(taxid)] <- 0
taxid[taxid == ''] <- 0
sample <- as.data.frame(total$V4)
sample <- sample %>% separate(`total$V4`, sep = '_', 
                              c('individual', 'time_point', 'sample'))

total <- cbind(total, taxonomy, sample, taxid)
colnames(total)[1:4] <- c('taxonomy', 'taxid', 'proportion', 'sample_info')
total[is.na(total)] <- 0
write.table(total, '~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/metaphlan/summarized_reports/20-4-1.mpa_summarized_total_report.txt', 
            col.names = T, row.names = F, quote = F, sep = '\t') # save the total report.

# phylum level stats
phylum <- total %>% filter(class == 0 & phylum != 0)
write.table(phylum$phylum_id, '../intermediates/20.taxonomy_profiling/metaphlan/mpa/phylum_id.txt', 
            sep = '\t', col.names = F, row.names = F, quote = F)
# taxonkit lineage phylum_id.txt | taxonkit reformat -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" > phylum_names.txt
phylum$phylum_new <- read.table('../intermediates/20.taxonomy_profiling/metaphlan/mpa/phylum_names.txt', 
                         sep = '\t', header = F, fill = NA)[, 4]
phylum[phylum == ''] <- 'unclassified'
phylum <- phylum %>% mutate(phylum_type = case_when(
  proportion >= 1 ~ phylum_new,
  proportion <1 & phylum_new != 'unclassified' ~ 'others',
  phylum_new == 'unclassified' ~ 'unclassified',
))
phylum$individual <- as.numeric(phylum$individual)
phylum$time_point <- as.numeric(phylum$time_point)
n_distinct(phylum$phylum_type)
color_8 <- c("#8BE380", "#DDD37A", "#E08883", "#A8D7CD", "#BC58D3", "#AA9FD7", "#bdbdbd", "#737373")
phylum$phylum_type <- factor(phylum$phylum_type, 
                             levels = c('p__Actinobacteria', 'p__Bacteroidetes', 
                                        'p__Euryarchaeota', 'p__Firmicutes', 
                                        'p__Proteobacteria', 'p__Verrucomicrobia',
                                        'others'))
phylum %>% ggplot(aes(x = time_point, y = proportion)) +
  geom_col(aes(fill = phylum_type)) + 
  facet_wrap(.~individual, nrow = 2) + 
  scale_fill_manual(values = color_8)+ 
  scale_x_continuous(limits = c(0.5, 10.5), 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) #~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/metaphlan/metaphlan_figures/20-4-1.phylum_level_composition
write.table(phylum, '~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/metaphlan/summarized_reports/20-4-2.mpa_phylum_report.txt', 
            col.names = T, row.names = F, quote = F, sep = '\t') # save the phylum-level summarized report.

# genus level stats
genus <- total %>% filter(genus != 0 & species == 0)
write.table(genus$genus_id, '../intermediates/20.taxonomy_profiling/metaphlan/mpa/genus_id.txt', 
            sep = '\t', col.names = F, row.names = F, quote = F)
# taxonkit lineage genus_id.txt | taxonkit reformat -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" > genus_names.txt
genus$genus_new <- read.table('../intermediates/20.taxonomy_profiling/metaphlan/mpa/genus_names.txt', 
                                sep = '\t', header = F, fill = NA)[, 8]
genus[genus == ''] <- 'unclassified'
genus <- genus %>% mutate(genus_type = case_when(
  proportion >= 2 ~ genus_new, 
  proportion < 2 & genus_new != 'unclassified' ~ 'others',
  genus_new == 'unclassified' ~ 'unclassified'
)) # the percentage of unclassified is much lower in mpa than kraken and mmseqs+bowtie, so the threshold of defining "others" is set 2% here, rather than 1%.
genus$individual <- as.numeric(genus$individual)
genus$time_point <- as.numeric(genus$time_point)
n_distinct(genus$genus_type)
color_32 <- c("#C8E8C1", "#A8A390", "#5ADC7F", "#DBB65D", "#E2888F", "#AECDE4", "#6A8EA8", "#E7C4A9",
              "#DDEAE4", "#DD47D9", "#70E0D9", "#B9E88C", "#D68A56", "#E55446", "#8A39EA", "#E54F9A",
              "#5BC0E8", "#A57099", "#E4DB55", "#E7E79E", "#C5F057", "#6389DC", "#83E6AF", "#74A16B",
              "#E99BD4", "#74E547", "#685BCB", "#B3A7E9", "#C774D8", "#E0C3D9", "#bdbdbd", "#737373")
genus %>% ggplot(aes(x = time_point, y = proportion)) + 
  geom_col(aes(fill = genus_type)) + 
  facet_wrap(.~ individual, nrow = 2) + 
  scale_fill_manual(values = color_32[4:32]) + 
  scale_x_continuous(limits = c(0.5, 10.5), 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/metaphlan/metaphlan_figures/20-4-2.genus_level_composition
write.table(genus, '~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/metaphlan/summarized_reports/20-4-3.mpa_genus_report.txt', 
            col.names = T, row.names = F, quote = F, sep = '\t') # save the genus-level summarized report.
