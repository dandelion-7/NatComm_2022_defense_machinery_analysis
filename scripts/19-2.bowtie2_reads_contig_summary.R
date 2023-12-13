# This script is for summarize the mapping results of CRISPR reads onto assembled contigs, so as to analyze the taxonomic composition of CRISPR loci recovered from metagenome.
# In the first part, mapping results that only reported one alignment even alignment time >1, were analyzed.
# In the second part, all alignments were reported into the bed file. In the analysis, first multi-aligned reads are evaluated for cross-contig alignment, if a read is repetitively aligned to one contig, then there will be no ambiguity in taxonomy assignment. Second, for reads with cross-contig align, if all these contigs belong to the same taxonomy, then there's no ambiguity in taxonomy assignment. Finally, if cross-aligned contigs span different taxonomy, if the dominant taxonomy is dramatically dominant (phylum level 100%, or genus level >90%), then this taxonomy is assigned to the contig. If the dominant taxonomy is not well prominent (phylum 50%-100%, genus 50%-90%), then 'uncertain' is assigned. If the dominant taxonomy is <50%, then unclassified is assigned.
#------------------------------------------------------------------------------------------------------------------------------------------

#install.packages('randomcoloR')
library(randomcoloR)
library(dplyr)
library(tidyr)
library(ggplot2)
getwd()
setwd("/home/zhanggaopu/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts")

# Preprocessing
bed <- read.table('../intermediates/19.bowtie2_reads_contig/bam/total.bed', 
                  sep = '\t', header = F)
reads <- as.data.frame(bed$V4 )
reads <- reads %>% separate(`bed$V4`, sep = '_', c('individual', 'time_point', 'sample', 'group', 'group_number', 'SRR', 'read'))
reads <- reads[, -4]
reads$read <- paste(reads$SRR, reads$read, sep = '_')
bed <- cbind(bed, reads)
colnames(bed)[1:6] <- c('contig', 'start', 'end', 'read_info', 'score', 'strand')
bed <- bed %>% mutate(aligned_length = end - start)
bed$contig_info <- paste(bed$individual, bed$contig, sep = '_')

bed %>% ggplot(aes(x = score)) + geom_histogram()
bed %>% ggplot(aes(x = score, y = aligned_length)) + 
  geom_point()

contigs <- read.table('../intermediates/15.mmseqs_taxonomy/taxonomy_new/mmseqs2_total_summary.txt', 
                      header=T, sep = '\t')
# Merging taxonomy annotation of contigs and alignments of reads to contigs.
contigs <- contigs[, 8:16]
bed <- bed[, c(5, 7, 8, 9, 10, 12, 14)]
contigs$contig_info <- paste(contigs$individual, contigs$contig_info, sep = '_')
total <- left_join(bed, contigs, by = "contig_info")
is.na(total) <- 'unidentified'

# Analyzing taxonomical composition of recovered CRISPRs.
crispr_genus <- total %>% group_by(individual.x, time_point, sample, group_number, genus) %>% summarise(read_count = n())
crispr_genus$group_info <- paste(crispr_genus$individual.x, crispr_genus$time_point, crispr_genus$sample, crispr_genus$group_number, sep = '_')
crispr_phylum <- total %>% group_by(individual.x, time_point, sample, group_number, phylum) %>% summarise(read_count = n())
crispr_phylum$group_info <- paste(crispr_phylum$individual.x, crispr_phylum$time_point, crispr_phylum$sample, crispr_phylum$group_number, sep = '_')

cross_genus_stat <- crispr_genus %>% filter(genus != 'unclassified' & genus != 'unidentified') %>% 
  group_by(individual.x, time_point, sample, group_number) %>% 
  summarise(genus_count = n(), read_count=sum(read_count))
cross_genus_stat %>% ggplot(aes(x = genus_count)) + geom_histogram() + 
  scale_x_continuous(limits = c(0, 10), breaks = c(0, 2, 4, 6, 8, 10)) + 
  theme_minimal() #19-2-3.cross_genus

cross_phylum_stat <- crispr_phylum %>% filter(phylum != 'unclassified' & phylum != 'unidentified') %>% 
  group_by(individual.x, time_point, sample, group_number) %>% 
  summarise(phylum_count = n(), read_count = sum(read_count))
cross_phylum_stat %>% ggplot(aes(x = phylum_count)) + geom_histogram() + 
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6)) + 
  theme_minimal() #19-2-4.cross_phylum

crispr_sum <- crispr_genus %>% group_by(individual.x, time_point, sample) %>% summarize(total_read_count = sum(read_count))

crispr_phylum <- left_join(crispr_phylum, crispr_sum)
crispr_phylum <- crispr_phylum %>% mutate(phylum_proportion = read_count/total_read_count)
crispr_phylum <- crispr_phylum %>% mutate(phylum_type = case_when(
  phylum == 'unclassified' ~ 'unclassified', 
  phylum_proportion >= 0.001 & phylum != 'unclassified' ~ phylum,
  phylum_proportion < 0.001 & phylum != 'unclassified' ~ 'others'
))



manual_color=c('#8dd3c7','#f16913','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#2171b5','#bc80bd','#ccebc5','#ffed6f',
               '#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#bdbdbd', '#7fc97f','#beaed4','#fdc086', '#d9d9d9', '#969696')
crispr_phylum$individual.x <- factor(crispr_phylum$individual.x, 
                                    levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
crispr_phylum$time_point <- factor(crispr_phylum$time_point, 
                                  levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
crispr_phylum %>% ggplot(aes(x = time_point, y = phylum_proportion)) + 
  facet_grid(. ~ individual.x) + 
  geom_col(aes(fill = phylum_type)) + 
  scale_fill_manual(values = distinctColorPalette(25)) + 
  theme(legend.position = 'top') + 
  guides(fill = guide_legend(nrow = 2)) + 
  theme(legend.text = element_text(size = 10)) #19-2-1.crispr_phylum

crispr_genus <- left_join(crispr_genus, crispr_sum)
crispr_genus <- crispr_genus %>% mutate(genus_proportion = read_count/total_read_count)
crispr_genus <- crispr_genus %>% mutate(genus_type = case_when(
  genus == 'unclassified' ~ 'unclassified',
  genus_proportion >= 0.01 & genus != 'unclassified' ~ genus,
  genus_proportion < 0.01  & genus != 'unclassified' ~ 'others'
))

manual_color_2 <- c('#c6dbef', '#fcbba1',	'#c7e9c0', '#dadaeb', '#2171b5', '#cb181d',	'#238b45', '#6a51a3',
                    '#9ecae1', '#fc9272', '#a1d99b', '#bcbddc', '#08519c', '#a50f15', '#006d2c', '#54278f',
                    '#6baed6', '#fb6a4a',	'#74c476', '#9e9ac8', '#08306b', '#67000d',	'#00441b', '#3f007d', 
                    '#4292c6', '#ef3b2c', '#41ab5d', '#807dba',
                    '#fbb4ae', '#b3cde3',	'#ccebc5', '#decbe4', '#fed9a6', '#ffffcc', '#e5d8bd', '#fddaec', '#f2f2f2',
                    '#8dd3c7', '#ffffb3',	'#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9',
                    '#e41a1c', '#377eb8',	'#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999')
crispr_genus$individual.x <- factor(crispr_genus$individual.x, 
                           levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
crispr_genus$time_point <- factor(crispr_genus$time_point, 
                                    levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
crispr_genus %>%  ggplot(aes(x = time_point, y = genus_proportion)) + 
  facet_grid(. ~ individual.x) + 
  geom_col(aes(fill = genus_type)) + 
  scale_fill_manual(values = distinctColorPalette(55)) + 
  theme(legend.position = 'top') + 
  guides(fill = guide_legend(nrow = 4)) + 
  theme(legend.text = element_text(size = 10)) #19-2-2.crispr_genus

# 19-3
# Compare the taxonomic composition of CRISPR reads and metagenomic reads.
# At the phylum level
bowtie2_phylum <- read.table("~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/mmseqs2+bowtie2/mmseqs2+bowtie2_summarized_reports/20-5-2.bowtie2_phylum_report.txt", 
                             header = T, sep = '\t')
bowtie2_phylum$object <- "total_reads"
bowtie2_phylum <- bowtie2_phylum[, c(1, 2, 4, 7, 8, 9)]

crispr_phylum$object <- "crispr_reads"
crispr_phylum <- crispr_phylum[, c(1, 2, 5, 9, 10, 11)]
colnames(crispr_phylum) <- colnames(bowtie2_phylum)
crispr_phylum$individual <- as.integer(crispr_phylum$individual)
crispr_phylum$time_point <- as.integer(crispr_phylum$time_point)

phylum <- rbind(crispr_phylum, bowtie2_phylum)
phylum_level <- c(sort(unique(phylum$phylum_type))[2:20], 'others', 'unclassified')
phylum$phylum_type <- factor(phylum$phylum_type, levels = phylum_level)
color_22 <- distinctColorPalette(19)
color_22[20:22] <- c("#f0f0f0", "#bdbdbd", "#737373") 
phylum %>% ggplot(aes(x = time_point, y = phylum_propotion)) + 
  geom_col(aes(fill = phylum_type)) + 
  facet_grid(object ~ individual) + 
  scale_fill_manual(values = color_22) +
  scale_x_continuous(limits = c(0.5, 10.5), 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
  guides(fill = guide_legend(ncol = 2, )) + 
  theme(legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size = 10))
#19-3-1.crispr&total_reads_phylum_comparison

# At the genus level
bowtie2_genus <- read.table("~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/mmseqs2+bowtie2/mmseqs2+bowtie2_summarized_reports/20-5-3.bowtie2_genus_report.txt", 
                            header = T, sep = '\t')
bowtie2_genus$object <- "total_reads"
bowtie2_genus <- bowtie2_genus[, c(1, 2, 4, 7, 8, 9)]

crispr_genus$object <- "crispr_reads"
crispr_genus <- crispr_genus[, c(1, 2, 5, 9, 10, 11)]
colnames(crispr_genus) <- colnames(bowtie2_genus)
crispr_genus$individual <- as.integer(crispr_genus$individual)
crispr_genus$time_point <- as.integer(crispr_genus$time_point)

genus <- rbind(crispr_genus, bowtie2_genus)
n_distinct(genus$genus_type)
color_65 <- distinctColorPalette(62)
color_65[63:65] <- c("#f0f0f0", "#bdbdbd", "#737373")
genus %>% ggplot(aes(x = time_point, y = genus_proportion)) + 
  geom_col(aes(fill = genus_type)) + 
  facet_grid(object ~ individual) + 
  scale_fill_manual(values = color_65) + 
  scale_x_continuous(limits = c(0.5, 10.5), 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
  guides(fill = guide_legend(ncol = 2, )) + 
  theme(legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size = 10))
# 19-3-2.crispr&total_reads_genus_comparison

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Part 2

# Re-analyze with bowtie2 reporting all the mapping results, and multiple alignment is summarized. 
bed <- read.table('../intermediates/19.bowtie2_reads_contig/bam_all/total.bed', 
                  sep = '\t', header = F)
reads <- as.data.frame(bed$V4)
reads <- reads %>% separate(`bed$V4`, sep = '_', c('individual', 'time_point', 'sample', 'group', 'group_number', 'SRR', 'read'))
reads <- reads[, -4]
reads$read <- paste(reads$SRR, reads$read, sep = '_')
bed <- cbind(bed, reads)
colnames(bed)[1:6] <- c('contig', 'start', 'end', 'read_info', 'score', 'strand')
bed <- bed %>% mutate(aligned_length = end - start)
bed$contig_info <- paste(bed$individual, bed$contig, sep = '_')

# Summarize multi-aligned reads (one read mapped to one contigs for multiple times)
repetitive_align_stat <- bed %>% group_by(contig_info, read_info) %>% 
  summarise(repetitve_align = n(), 
            max_length = max(aligned_length), 
            min_length = min(aligned_length))
repetitive_align_stat %>% ggplot(aes(x = repetitve_align)) + 
  geom_histogram(binwidth = 0.5) +
  scale_x_continuous(limits = c(0.5, 55))

contigs <- read.table('../intermediates/15.mmseqs_taxonomy/taxonomy_new/mmseqs2_total_summary.txt', 
                      header=T, sep = '\t')
# Merging taxonomic annotation of contigs with alignments of reads to contigs.
contigs <- contigs[, 8:16]
contigs$contig_info <- paste(contigs$individual, contigs$contig_info, sep = '_')
contigs <- contigs[, -1]
total <- left_join(repetitive_align_stat, contigs, by = "contig_info")

getmode <- function(x){
  uniq <- unique(x)
  freq <- tabulate(match(x, uniq))
  mf <- max(freq)
  uniq[freq == mf]
}
getmodefreq <- function(x){
  uniq <- unique(x)
  freq <- tabulate(match(x, uniq))
  mf <- max(freq)
  mf/length(x)
}
getmode(c('a', 'b', 'c', 'c', 'c', 'd', 'd', 'd'))[1]
getmodefreq(c('a', 'b', 'c', 'c', 'c', 'd', 'd', 'd'))


cross_contig_align_stat <- total %>% group_by(read_info) %>% 
  summarise(cross_aligned_contig = n(), 
            cross_aligned_phylum = n_distinct(phylum), 
            cross_aligned_genus = n_distinct(genus),
            dominant_phylum = getmode(phylum)[1],
            dominant_phylum_freq = getmodefreq(phylum),
            dominant_genus = getmode(genus)[1],
            dominant_genus_freq = getmodefreq(genus),
            total_align = sum(repetitve_align))
cross_contig_align_stat %>% ggplot(aes(x = cross_aligned_contig)) + 
  geom_histogram(binwidth = 0.5) + 
  scale_x_continuous(limits = c(0.5, 40)) # most reads are only mapped to one contig.
cross_contig_align_stat %>% filter(cross_aligned_contig > 1) %>% ggplot(aes(x = cross_aligned_genus)) +
  geom_histogram()
cross_contig_align_stat %>% filter(cross_aligned_contig > 1) %>%  ggplot(aes(x = cross_aligned_phylum)) +
  geom_histogram()

cross_contig_align_stat %>% filter(cross_aligned_contig >= 1) %>%  ggplot(aes(x = dominant_phylum_freq)) +
  geom_histogram()
cross_contig_align_stat %>% filter(cross_aligned_contig > 1) %>% ggplot(aes(x = dominant_genus_freq)) +
  geom_histogram()

total_reads_stat <- cross_contig_align_stat %>% mutate(
  phylum = case_when(
  dominant_phylum == 'unclassified' ~ 'unclassified',
  dominant_phylum_freq == 1 ~ dominant_phylum,
  dominant_phylum_freq < 1 & dominant_phylum_freq >= 0.5 & dominant_phylum != 'unclassified' ~ 'uncertain',
  dominant_phylum_freq < 0.5 ~ 'unclassified'), 
  genus = case_when(
  dominant_genus == 'unclassified' ~ 'unclassified',
  dominant_genus_freq >= 0.9 ~ dominant_genus,
  dominant_genus_freq < 0.9 & dominant_genus_freq >= 0.5 & dominant_genus != 'unclassified' ~ 'uncertain',
  dominant_genus_freq < 0.5 ~ 'unclassified')
)
# For reads cross-aligned to contigs spanning different taxonomy, according to the abundancy of dominant taxonomy, the reads are assigned with the taxonomy.

reads <- as.data.frame(total_reads_stat$read_info)
reads <- reads %>% separate(`total_reads_stat$read_info`, sep = '_', c('individual', 'time_point', 'sample', 'group', 'group_number', 'SRR', 'read'))
total_reads_stat <- cbind(total_reads_stat, reads[, c(1, 2, 3, 5)])
total_reads_stat$sample_info <- paste(total_reads_stat$individual, total_reads_stat$time_point, total_reads_stat$sample, sep = '_')
sample_read_sum <- total_reads_stat %>% group_by(sample_info) %>% 
  summarise(read_count = n())

total_reads_stat %>% write.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/19.bowtie2_reads_contig/figures/19-3-1.crispr_reads_taxonomy_stat.txt',
                                 col.names = T, row.names = F, quote = F, sep = '\t')
total_reads_stat <- read.delim('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/19.bowtie2_reads_contig/figures/19-3-1.crispr_reads_taxonomy_stat.txt',
                               header = T, sep = '\t')

# At the phylum-level, compare taxonomic compositions of crispr reads and total metagenomic reads.
total_reads_phylum_stat <- total_reads_stat %>% group_by(individual, time_point, sample, sample_info, phylum) %>% 
  summarise(phylum_count = n())
total_reads_phylum_stat <- left_join(total_reads_phylum_stat, sample_read_sum, by = 'sample_info')
total_reads_phylum_stat <- total_reads_phylum_stat %>% mutate(
  phylum_proportion = phylum_count/read_count, 
  phylum_type = case_when(
    !(str_detect(phylum, 'p__')) ~ phylum,
    phylum_count/read_count < 0.01 ~ 'others',
    phylum_count/read_count >= 0.01 ~ phylum
  ))
total_reads_phylum_stat$object <- 'CRISPR_reads'
total_reads_phylum_stat <- total_reads_phylum_stat[, c(1, 2, 4, 5, 8, 9, 10)]

bowtie2_phylum <- read.table("~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/mmseqs2+bowtie2/mmseqs2+bowtie2_summarized_reports/20-5-2.bowtie2_phylum_report.txt", 
                             header = T, sep = '\t')
bowtie2_phylum$object <- "microbiome_reads"
bowtie2_phylum <- bowtie2_phylum[, c(1, 2, 3, 4, 7, 8, 9)]
colnames(total_reads_phylum_stat) <- colnames(bowtie2_phylum)
total_reads_phylum_stat$individual <- as.integer(total_reads_phylum_stat$individual)
total_reads_phylum_stat$time_point <- as.integer(total_reads_phylum_stat$time_point)
phylum_comparison <- rbind(bowtie2_phylum, total_reads_phylum_stat)
n_distinct(phylum_comparison$phylum_type)
sort(unique(phylum_comparison$phylum_type))
phylum_comparison$phylum_type <- factor(phylum_comparison$phylum_type, 
                                        levels = c(sort(unique(phylum_comparison$phylum_type))[2:10], 'others', 'uncertain', 'unclassified'))
color_13 <- distinctColorPalette(9)
color_13[10:13] <- c('#ffffff', '#d9d9d9', '#969696', '#737373')
phylum_comparison %>% ggplot(aes(x = time_point, y = phylum_propotion)) + 
  geom_col(aes(fill = phylum_type)) + 
  facet_grid(object ~ individual) + 
  scale_fill_manual(values = color_13) +
  scale_x_continuous(limits = c(0.5, 10.5), 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
  guides(fill = guide_legend(ncol = 1, )) + 
  theme(legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size = 10))
# ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/19.bowtie2_reads_contig/figures/19-3-1.crispr_total_reads_phylum_comparison_multialign_summarized

# At the genus-level, compare taxonomic compositions of crispr reads and total metagenomic reads.
total_reads_genus_stat <- total_reads_stat %>% group_by(individual, time_point, sample, sample_info, genus) %>% 
  summarise(genus_count = n())
total_reads_genus_stat <- left_join(total_reads_genus_stat, sample_read_sum, by = 'sample_info')
total_reads_genus_stat <- total_reads_genus_stat %>% mutate(
  genus_proportion = genus_count/read_count, 
  genus_type = case_when(
    !(str_detect(genus, 'g__')) ~ genus,
    genus_count/read_count < 0.01 & str_detect(genus, 'g__') ~ 'others',
    genus_count/read_count >= 0.01 ~ genus
  ))
total_reads_genus_stat$object <- 'CRISPR_reads'
total_reads_genus_stat <- total_reads_genus_stat[, c(1, 2, 4, 5, 8, 9, 10)]

bowtie2_genus <- read.table("~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/mmseqs2+bowtie2/mmseqs2+bowtie2_summarized_reports/20-5-3.bowtie2_genus_report.txt", 
                             header = T, sep = '\t')
bowtie2_genus$object <- "microbiome_reads"
bowtie2_genus <- bowtie2_genus[, c(1, 2, 3, 4, 7, 8, 9)]
colnames(total_reads_genus_stat) <- colnames(bowtie2_genus)
total_reads_genus_stat$individual <- as.integer(total_reads_genus_stat$individual)
total_reads_genus_stat$time_point <- as.integer(total_reads_genus_stat$time_point)
genus_comparison <- rbind(total_reads_genus_stat, bowtie2_genus)
n_distinct(genus_comparison$genus_type)
color_63 <- distinctColorPalette(61)
color_63[60:63] <- c('#ffffff', '#d9d9d9', '#969696', '#737373')
genus_comparison %>% ggplot(aes(x = time_point, y = genus_proportion)) + 
  geom_col(aes(fill = genus_type)) + 
  facet_grid(object ~ individual) + 
  scale_fill_manual(values = color_63) + 
  scale_x_continuous(limits = c(0.5, 10.5), 
                     breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
  guides(fill = guide_legend(ncol = 2, )) + 
  theme(legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size = 10))
# ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/19.bowtie2_reads_contig/figures/19-2-2.crispr_total_reads_genus_comparison_multialign_summarized


# Summarize cross-phylum/genus distribution of CRISPR groups, not finished.
total_reads_stat$group_info <- paste(total_reads_stat$sample_info, total_reads_stat$group_number, sep = '_')
group_read_sum <- total_reads_stat %>% group_by(group_info) %>% 
  summarize(group_read_count = n())
total_reads_stat <- left_join(total_reads_stat, group_read_sum, by = 'group_info')
total_reads_stat %>% write.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/19.bowtie2_reads_contig/figures/19-3-1.crispr_reads_taxonomy_stat.txt',
                                 col.names = T, row.names = F, quote = F, sep = '\t')

large_groups <- total_reads_stat %>% filter(group_read_count >= 5)
large_groups_cross_genus_stats <- large_groups %>% filter(genus != 'unclassified' & genus != 'uncertain')
classified_reads_genus_stat <- large_groups_cross_genus_stats %>% group_by(group_info) %>% 
  summarize(classified_read_count = n())
large_groups_cross_genus_stats <- left_join(large_groups_cross_genus_stats, classified_reads_genus_stat)
large_groups_cross_genus_stats <- large_groups_cross_genus_stats %>% mutate(
  classified_proportion = classified_read_count/group_read_count
)
cross_genus_stat <- large_groups_cross_genus_stats %>% group_by(group_info, group_read_count, classified_read_count, classified_proportion) %>% 
  summarise(spanned_genus = n_distinct(genus), 
            dominant_genus = getmode(genus)[1], 
            dominant_genus_freq = getmodefreq(genus), 
            sparse_genus = getsparse(genus)[1],
            sparse_genus_freq = getsparsefreq(genus)
            )
cross_genus_stat %>% ggplot(aes(x = spanned_genus)) + 
  geom_histogram(binwidth = 0.5)

getsparse <- function(x){
  uniq <- unique(x)
  freq <- tabulate(match(x, uniq))
  mf <- min(freq)
  uniq[freq == mf]
}
getsparse(c('a', 'a', 'b', 'c', 'd', 'd', 'd'))[1]

getsparsefreq <- function(x){
  uniq <- unique(x)
  freq <- tabulate(match(x, uniq))
  mf <- min(freq)
  mf/length(x)
}
getsparsefreq(c('a', 'b', 'c', 'c', 'c', 'd', 'd', 'd'))
