# This script is for summarizing the taxonomic annotation results from Kraken2, following script 20-1.
#----------------------------------------------------------------------------------------------------------------------
library(dplyr)
library(randomcoloR)
getwd()
setwd("~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts")

# pre-processing
total <- read.table('../intermediates/20.taxonomy_profiling/kraken/total_report.txt', 
                    header = F, sep = '\t')
n_distinct(total$V7)
sample_info <- as.data.frame(total$V7)
sample_info <- sample_info %>% separate(`total$V7`, sep = '_', into = c('individual', 'time_point', 'sample'))
total <- cbind(total, sample_info)
colnames(total)[1:7] <- c('proportion', 'included_reads', 'current_reads', 'level', 'taxid', 'taxonomy', 'sample_info')

total_reads <- total %>% group_by(sample_info) %>% 
  summarise(total_reads = sum(current_reads))
total <- left_join(total, total_reads)
total <- total %>% mutate(included_proportion = included_reads / total_reads, 
                          current_proportion = current_reads / total_reads)
total <- total %>% mutate(level_number = case_when(
  str_detect(level, 'U') ~ 0,
  str_detect(level, 'R') ~ 1,
  str_detect(level, 'D') ~ 2,
  str_detect(level, 'P') ~ 3,
  str_detect(level, 'C') ~ 4,
  str_detect(level, 'O') ~ 5,
  str_detect(level, 'F') ~ 6,
  str_detect(level, 'G') ~ 7,
  str_detect(level, 'S') ~ 8
))
write.table(total, '~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/kraken/summarized_reports/20-3-1.kraken_summarized_total_report.txt', 
            sep = '\t', col.names = T, quote = F, row.names = F)

# phylum-level stats
phylum <- total %>% filter(level_number < 3 | level == 'P')
phylum <- phylum %>% mutate(phylum_taxonomy = case_when(
  level == 'P' ~ paste('p__', taxonomy, sep = ''),
  level != 'P' ~ 'unclassified'
), phylum_proportion = case_when(
  level == 'P' ~ included_reads/total_reads, 
  level != 'P' ~ current_reads/total_reads
))
phylum <- phylum %>% mutate(phylum_type = case_when(
  phylum_taxonomy == 'unclassified' ~ 'unclassified',
  phylum_proportion >= 0.01 ~ phylum_taxonomy, 
  phylum_proportion < 0.01 ~ 'others'
))
phylum$individual <- as.numeric(phylum$individual)
phylum$time_point <- as.numeric(phylum$time_point)
n_distinct(phylum$phylum_type)
#color_8 <- distinctColorPalette(6)
#color_8[7:8] <- c("#bdbdbd", "#737373")
color_8 <- c("#8BE380", "#DDD37A", "#E08883", "#A8D7CD", "#BC58D3", "#AA9FD7", "#bdbdbd", "#737373")
phylum$phylum_type <- factor(phylum$phylum_type, 
                             levels = c('p__Actinomycetota', 'p__Bacillota', 
                                        'p__Bacteroidota', 'p__Euryarchaeota', 
                                        'p__Pseudomonadota', 'p__Verrucomicrobiota', 
                                        'others', 'unclassified'))
phylum %>% ggplot(aes(x = time_point, y = phylum_proportion)) + 
  geom_col(aes(fill = phylum_type)) + 
  facet_wrap(.~individual, nrow = 2) + 
  scale_fill_manual(values = color_8) + 
  scale_x_continuous(limits = c(0.5, 10.5), 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) # ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/kraken/kraken_figures/20-3-1.phylum_level_composition
write.table(phylum, '~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/kraken/summarized_reports/20-3-2.kraken_phylum_reports.txt', 
            sep = '\t', col.names = T, row.names = F, quote = F) # save the phylum-level table

# genus-level stats
genus <- total %>% filter(level_number < 7 | level == 'G')
genus <- genus %>% mutate(genus_taxonomy = case_when(
  level == 'G' ~ paste('g__', taxonomy, sep = ''),
  level != 'G' ~ 'unclassified'
), genus_proportion = case_when(
  level == 'G' ~ included_reads/total_reads, 
  level != 'G' ~ current_reads/total_reads
))

genus <- genus %>% mutate(genus_type = case_when(
  genus_taxonomy == 'unclassified' ~ 'unclassified',
  genus_proportion >= 0.01 ~ genus_taxonomy, 
  genus_proportion < 0.01 ~ 'others'
))
genus$individual <- as.numeric(genus$individual)
genus$time_point <- as.numeric(genus$time_point)
n_distinct(genus$genus_type)
#color_32 <- distinctColorPalette(30)
#color_32[31:32] <- c("#bdbdbd", "#737373")
color_32 <- c("#C8E8C1", "#A8A390", "#5ADC7F", "#DBB65D", "#E2888F", "#AECDE4", "#6A8EA8", "#E7C4A9",
              "#DDEAE4", "#DD47D9", "#70E0D9", "#B9E88C", "#D68A56", "#E55446", "#8A39EA", "#E54F9A",
              "#5BC0E8", "#A57099", "#E4DB55", "#E7E79E", "#C5F057", "#6389DC", "#83E6AF", "#74A16B",
              "#E99BD4", "#74E547", "#685BCB", "#B3A7E9", "#C774D8", "#E0C3D9", "#bdbdbd", "#737373")
genus %>% ggplot(aes(x = time_point, y = genus_proportion)) + 
  geom_col(aes(fill = genus_type)) + 
  facet_wrap(.~individual, nrow = 2) + 
  scale_fill_manual(values = color_32) + 
  scale_x_continuous(limits = c(0.5, 10.5), 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) # ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/kraken/kraken_figures/20-3-2.genus_level_composition
genus %>% filter(genus_type != 'unclassified') %>% ggplot(aes(x = time_point, y = genus_proportion)) + 
  geom_col(aes(fill = genus_type), color = '#d9d9d9') + 
  facet_wrap(.~individual, nrow = 2) + 
  scale_fill_manual(values = color_32) + 
  scale_x_continuous(limits = c(0.5, 10.5), 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) # ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/kraken/kraken_figures/20-3-3.genus_level_composition_no_unclassified
write.table(genus, '~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/kraken/summarized_reports/20-3-3.kraken_genus_reports.txt', 
            sep = '\t', col.names = T, row.names = F, quote = F)
