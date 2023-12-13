# abundance analysis of defense systems are normalized to gene length and system size in this script.

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(randomcoloR)

getwd()
setwd("/home/zhanggaopu/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts")

# ==============================================================================
# summarize the total mapping count table
total_counts <- read.delim('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/featureCounts_new/merged_gene_count.txt', 
                           sep = '\t', header = F)
colnames(total_counts) <- c('gene_name', 'CDS', '10', '1', '2', '3', '4', '5', 
                            '6', '7', '8', '9', 'subject')
total_counts <- total_counts %>% gather(key = 'time_point', value = 'mapping_count', 
                                        colnames(total_counts)[3:12])
total_counts$time_point <- as.numeric(total_counts$time_point)
sample_read_sum <- total_counts %>% group_by(subject, time_point) %>% 
  summarise(sample_read_count = sum(mapping_count))
total_counts <- left_join(total_counts, sample_read_sum, by = c('subject', 'time_point'))
defense_gene_counts <- total_counts %>% filter(gene_name != 'non-defense-gene')

total_counts %>% write.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/defense_gene_stats/26-4-1.summarized_total_counts.txt', 
                             sep = '\t', row.names = F, col.names = T, quote = F) # write table
defense_gene_counts %>% write.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/defense_gene_stats/26-4-2.summarized_defense_gene_counts.txt',
                                    sep = '\t', row.names = F, col.names = T, quote = F) # write table
sample_read_sum %>% write.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/defense_gene_stats/26-4-3.sample_read_sum.txt', 
                                sep = '\t', row.names = F, col.names = T)
# ==============================================================================

# ==============================================================================
# summarize the mapping count of predicted defense genes
defense_gene_counts <- read.delim('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/defense_gene_stats/26-4-2.summarized_defense_gene_counts.txt', 
                                  sep = '\t', header = T)
sample_defense_gene_read_sum <- defense_gene_counts %>% group_by(subject, time_point) %>% 
  summarise(sample_defense_gene_read_count = sum(mapping_count))
defense_gene_counts <- left_join(defense_gene_counts, sample_defense_gene_read_sum, by = c('subject', 'time_point'))
# ==============================================================================

# ==============================================================================
#summarize the proportion of defense gene read counts among the total metagenomic reads.
sample_read_sum <- read.delim('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/defense_gene_stats/26-4-3.sample_read_sum.txt', 
                              sep = '\t', header = T)
gene_counts <- left_join(sample_read_sum, sample_defense_gene_read_sum, by = c('subject', 'time_point'))
gene_counts %>% ggplot(aes(x = time_point)) + 
  geom_col(aes(y = sample_defense_gene_read_count / sample_read_count * 100), width = 0.8, fill = '#bdbdbd') +
  geom_line(aes(y = sample_read_count / 1e+8)) +
  geom_point(aes(y = sample_read_count / 1e+8)) + 
  facet_wrap(subject~., nrow = 2) + 
  scale_y_continuous(name = "Predicted defense gene reads proportion(%)", sec.axis = sec_axis(~.*1e+8, name = 'Total reads')) + 
  scale_x_continuous(name = "Time point", breaks = c(1:10)) + 
  theme(axis.text.y.left = element_text(face = 'bold', color = '#bdbdbd', size = 12),
        axis.ticks.y.left = element_line(color = '#bdbdbd'),
        axis.text.y.right = element_text(face = 'bold', color = '#000000'), 
        axis.title.y.left = element_text(color = '#bdbdbd', face = 'bold'))
#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures_new/26-4-1.defense_gene_proportions
gene_counts %>% ggplot(aes(x = sample_read_count, 
                           y = sample_defense_gene_read_count/sample_read_count * 100)) + 
  geom_point() + 
  # geom_smooth(method = 'lm') + 
  labs(x = 'Total reads', y = "Predicted defense gene reads proportion(%)")
summary(lm(formula = gene_counts$sample_defense_gene_read_count ~ gene_counts$sample_read_count))
#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures_new/26-4-2.defense_gene_proportions_vs_total_reads
# ==============================================================================

# ==============================================================================
# read the annotation table of all predicted defense genes and systems.
defense_gene_annotation <- read.delim('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/gene_annotation/26-1-4.total_defense_gene_with_fine_annotation.txt', 
                                      sep = '\t', header = T)
defense_gene_annotation$gene_len <- as.numeric(str_remove_all(defense_gene_annotation$gene_len, pattern = "gene_len="))
colnames(defense_gene_annotation)
defense_gene_annotation <- defense_gene_annotation[, -c(31, 40, 41, 43, 44, 45)] #remove some columns not required for following analysis

# merge the defense gene annotation and count tables.
defense_genes <- left_join(defense_gene_counts, defense_gene_annotation, by = c("subject", 'CDS', 'gene_name'))
# ==============================================================================

# ==============================================================================
# summarize the average number of genes within a system of each type.
gene_number_in_type <- defense_genes %>% group_by(subject, time_point, type) %>% 
  summarize(gene_number = n()) # summarize the total number of genes belonging to one type of system, in each sample (subject+time_point)
system_number_type_level <- defense_genes %>% group_by(subject, time_point, type, sys_id.x) %>% 
  summarize(system_number = 1)
system_number_type_level <- system_number_type_level %>% group_by(subject, time_point, type) %>% 
  summarize(system_number = sum(system_number)) # summarize the number of systems of different types, in each sample (subject+time_point)
average_gene_number_in_type <- left_join(gene_number_in_type, system_number_type_level, by = c('subject', 'time_point', 'type'))
average_gene_number_in_type$average_system_size <- average_gene_number_in_type$gene_number / average_gene_number_in_type$system_number
# ==============================================================================

# ==============================================================================
# compute the different classes of proportions of all predicted defense genes.
defense_genes <- left_join(defense_genes, average_gene_number_in_type[, c(1, 2, 3, 6)],
                           by = c('subject', 'time_point', 'type'))
colnames(defense_genes)
defense_genes <- defense_genes %>% 
  mutate(absolute_proportion_of_gene = mapping_count/sample_read_count)
defense_genes <- defense_genes %>% 
  mutate(relative_proportion_of_gene = mapping_count/sample_defense_gene_read_count)
defense_genes <- defense_genes %>% 
  mutate(normalized_absolute_proportion_of_gene = mapping_count/sample_read_count/gene_len/average_system_size)
defense_genes <- defense_genes %>% 
  mutate(normalized_relative_proportion_of_gene = mapping_count/sample_defense_gene_read_count/gene_len/average_system_size)
# Four types of proportions are computed for each gene.
# ==============================================================================

# ==============================================================================
# Sum the proportions of genes at the type level.
defense_types <- defense_genes %>% group_by(subject, time_point, type) %>% 
  summarize(absolute_proportion_of_type = sum(absolute_proportion_of_gene), 
            relative_proportion_of_type = sum(relative_proportion_of_gene), 
            normalized_absolute_proportion_of_type = sum(normalized_absolute_proportion_of_gene), 
            normalized_relative_proportion_of_type = sum(normalized_relative_proportion_of_gene))
defense_types_sum <- defense_types %>% group_by(subject, time_point) %>% 
  summarize(absolute_proportion_sum = sum(absolute_proportion_of_type), 
            relative_proportion_sum = sum(relative_proportion_of_type), 
            normalized_absolute_proportion_sum = sum(normalized_absolute_proportion_of_type), 
            normalized_relative_proportion_sum = sum(normalized_relative_proportion_of_type))
# ==============================================================================

# ==============================================================================
# Merge the sum of all types' proportions to the list containing each types proportions
defense_types <- left_join(defense_types, defense_types_sum, by = c("subject", "time_point"))
# ==============================================================================

# ==============================================================================
# Use normalized relative proportions to calculate the abundances of defense system types.
defense_types_norm_relative <- defense_types[, c(1, 2, 3, 7, 11)]
defense_types_norm_relative_spread <- defense_types_norm_relative %>% spread(key = type, value = normalized_relative_proportion_of_type, fill = 0)
defense_types_norm_relative <- defense_types_norm_relative_spread %>% 
  gather(key = 'type', value = 'normalized_relative_proportion_of_type', colnames(defense_types_norm_relative_spread)[4:123]) 
# Through spread and gather, fill 0 if a specific type isn't present in a sample (subject+time_point).
defense_types_norm_relative <- defense_types_norm_relative %>% 
  mutate(type_abundance = normalized_relative_proportion_of_type/normalized_relative_proportion_sum)
average_norm_relative <- defense_types_norm_relative %>% group_by(type) %>% 
  summarize(average_type_abundance = sum(type_abundance)/100) # calculate the average abundances of different types among all the samples.
defense_types_norm_relative <- left_join(defense_types_norm_relative, average_norm_relative, by = 'type')
defense_types_norm_relative <- defense_types_norm_relative %>% 
  mutate(legend_average = case_when(
    average_type_abundance >= 0.005 ~ type,
    average_type_abundance < 0.005 ~ 'others'),
    legend_in_sample = case_when(
      type_abundance >= 0.02 ~ type,
      type_abundance < 0.02 ~ 'others'
    ))
# ==============================================================================

# ==============================================================================
# show the (average) relative abundances of different types of systems.
defense_types_norm_relative %>% 
  filter(average_type_abundance >= 0.01) %>% 
  ggplot(aes(x = reorder(type, -average_type_abundance), y = type_abundance * 100)) +
  geom_boxplot(width = 0.8, outlier.alpha = 0) +
  stat_boxplot(geom = 'errorbar', width = 0.3) +
  geom_jitter(aes(fill = as.character.numeric_version(subject)), shape = 21, width = 0.15, alpha = 0.6) +
  stat_summary(fun = 'mean', geom = 'point', shape = 25, color = 'red', fill = 'red', size = 2) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 10), 
        axis.text.y = element_text(size = 10)) + 
  scale_y_continuous(breaks = c(0, 1, 5, 10, 15, 20, 25, 30, 35)) + 
  theme(legend.position = 'None') +
  labs(x = 'Defense systems (average abundance > 1%)', 
       y = 'Designated defense system abundance in samples (%)') + 
  geom_hline(aes(yintercept = 1), linetype = 'dashed', color = '#737373')
#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures_new/26-4-3.system_average_abundances
# ==============================================================================

# ==============================================================================
# show the composition of defense systems in each sample.
n_distinct(defense_types_norm_relative$legend)
color_31 <- c('#bdbdbd', distinctColorPalette(30))
n_distinct(defense_types_norm_relative$legend_in_sample)
color_43 <- c('#bdbdbd', distinctColorPalette(42))
defense_types_norm_relative %>% ggplot(aes(x = as.character.numeric_version(time_point), 
                                           y = type_abundance)) + 
  geom_col(aes(fill = reorder(legend_in_sample, type_abundance))) + 
  facet_wrap(subject ~., nrow = 2) + 
  scale_fill_manual(values = color_43) + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = 'Time point', y = 'Relative abundance', fill = 'Defense system type') + 
  theme(legend.key.size = unit(4, 'mm'))
#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures_new/26-4-4.type_abundances_legend_sample
defense_types_norm_relative %>% ggplot(aes(x = as.character.numeric_version(time_point), 
                                           y = type_abundance)) + 
  geom_col(aes(fill = reorder(legend_average, average_type_abundance))) + 
  facet_wrap(subject ~., nrow = 2) + 
  scale_fill_manual(values = color_43) + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = 'Time point', y = 'Relative abundance', fill = 'Defense system type') + 
  theme(legend.key.size = unit(4, 'mm'))
#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures_new/26-4-5.type_abundances_legend_average
defense_types_norm_relative %>% filter(type != 'RM' & type != 'CasFinder') %>% 
  ggplot(aes(x = as.character.numeric_version(time_point), 
                                           y = type_abundance)) + 
  geom_col(aes(fill = reorder(legend_in_sample, type_abundance))) + 
  facet_wrap(subject ~., nrow = 2) + 
  scale_fill_manual(values = color_43) + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = 'Time point', y = 'Relative abundance', fill = 'Defense system type') + 
  theme(legend.key.size = unit(4, 'mm'))
#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures_new/26-4-6.type_abundances_legend_sample_no_RM_Cas
defense_types_norm_relative %>% filter(type != 'RM' & type != 'CasFinder') %>% 
  ggplot(aes(x = as.character.numeric_version(time_point), 
                                           y = type_abundance)) + 
  geom_col(aes(fill = reorder(legend_average, average_type_abundance))) + 
  facet_wrap(subject ~., nrow = 2) + 
  scale_fill_manual(values = color_43) + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(x = 'Time point', y = 'Relative abundance', fill = 'Defense system type') + 
  theme(legend.key.size = unit(4, 'mm'))
#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures_new/26-4-7.type_abundances_legend_average_no_RM_Cas
# ==============================================================================

# ==============================================================================
# combine the taxonomic annotation of contigs with total defense genes' information
total_contig_taxon <- read.delim('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/15.mmseqs_taxonomy/taxonomy_new/mmseqs2_total_summary.txt', 
                          sep = '\t', header = T)
total_contig_id <- as.data.frame(total_contig_taxon$contig)
total_contig_id <- separate(total_contig_id, col = `total_contig_taxon$contig`, sep = "_", 
                            into = c('subject', 'kmer', 'id'))
total_contig_id$contig_id <- paste(total_contig_id$kmer, total_contig_id$id, sep = '_')
total_contig_taxon <- cbind(total_contig_taxon, total_contig_id[, c(1, 4)])
total_contig_taxon$subject <- as.numeric(total_contig_taxon$subject)
defense_genes <- left_join(defense_genes, total_contig_taxon, by = c("subject", "contig_id"))
defense_genes %>% write.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/defense_gene_stats/26-4-4.defense_gene_counts_with_annotation&taxonomy.txt',
                                    sep = '\t', row.names = F, col.names = T, quote = F) # write table
# ==============================================================================

# ==============================================================================
# summarize the relative abundances of each cas gene with the normalized relative abundances
cas_genes <- defense_genes %>% filter(type == 'CasFinder')
cas_gene_names <- as.data.frame(cas_genes$gene_name)
cas_gene_names <- separate(cas_gene_names, col = `cas_genes$gene_name`, sep = "_", into = c('gene_name_simplified'))
n_distinct(cas_gene_names$gene_name_simplified)
cas_genes$gene_name_simplified <- cas_gene_names$gene_name_simplified
cas_genes_proportions <- cas_genes %>% group_by(subject, time_point, gene_name_simplified) %>% 
  summarize(absolute_proportion_sum = sum(absolute_proportion_of_gene), 
            relative_proportion_sum = sum(relative_proportion_of_gene), 
            normalized_absolute_proportion_of_gene = sum(normalized_absolute_proportion_of_gene), 
            normalized_relative_proportion_of_gene = sum(normalized_relative_proportion_of_gene))
cas_genes_norm_relative <- cas_genes_proportions[, c(1, 2, 3, 7)]
cas_genes_norm_relative_matrix <- cas_genes_norm_relative %>% spread(key = gene_name_simplified, value = normalized_relative_proportion_of_gene, fill = 0)
cas_genes_norm_relative <- cas_genes_norm_relative_matrix %>% 
  gather(key = 'gene_name_simplified', value = 'normalized_relative_proportion_of_cas_gene', 
         colnames(cas_genes_norm_relative_matrix)[3:59])
cas_genes_sum <- cas_genes_norm_relative %>% group_by(subject, time_point) %>% 
  summarize(normalized_cas_gene_relative_proportion_sum = sum(normalized_relative_proportion_of_cas_gene))
cas_genes_norm_relative <- left_join(cas_genes_norm_relative, cas_genes_sum, by = c("subject", "time_point"))
cas_genes_norm_relative <- cas_genes_norm_relative %>% 
  mutate(cas_gene_abundance = normalized_relative_proportion_of_cas_gene/normalized_cas_gene_relative_proportion_sum)
cas_genes_average_abundances <- cas_genes_norm_relative %>% group_by(gene_name_simplified) %>% 
  summarize(cas_gene_average_abundance = sum(cas_gene_abundance)/100)
cas_genes_norm_relative <- left_join(cas_genes_norm_relative, cas_genes_average_abundances, by = "gene_name_simplified")
cas_genes_norm_relative %>% filter(cas_gene_average_abundance >= 0.01) %>% 
  ggplot(aes(x = reorder(gene_name_simplified, -cas_gene_average_abundance), y = cas_gene_abundance)) + 
  geom_boxplot(width = 0.8, outlier.alpha = 0) +
  stat_boxplot(geom = 'errorbar', width = 0.3) +
  geom_jitter(aes(fill = as.character.numeric_version(subject)), shape = 21, width = 0.15, alpha = 0.6) + 
  scale_fill_manual(values = color_31)
# ==============================================================================

# ==============================================================================
# summarize the taxonomic distributions (phylum) of different types of defense systems, with normalized absolute proportion of the genes (RPKM)
type_phylum <- defense_genes %>% group_by(type, phylum) %>% 
  summarize(type_norm_abso_proportion_sum = sum(normalized_absolute_proportion_of_gene))
type_phylum <- drop_na(type_phylum)
type_phylum$type_norm_abso_proportion_sum <- type_phylum$type_norm_abso_proportion_sum * 1e9
type_phylum_matrix <- type_phylum %>% spread(key = phylum, value = type_norm_abso_proportion_sum, fill = 0)
type_phylum <- type_phylum_matrix %>% gather(key = 'phylum', value = 'type_norm_abso_proportion_sum', colnames(type_phylum_matrix)[2:19])
type_phylum_sum <- type_phylum %>% group_by(phylum) %>% summarize(phylum_norm_abso_proportion_sum = sum(type_norm_abso_proportion_sum))
type_phylum <- left_join(type_phylum, type_phylum_sum, by = 'phylum')
main_phylum <- c("p__Actinomycetota", "p__Bacillota", "p__Bacteroidota", "p__Pseudomonadota", "unclassified")
main_types <- unique(defense_types_norm_relative[defense_types_norm_relative$average_type_abundance >= 0.01,][['type']])
type_phylum %>% filter(phylum %in% main_phylum & type %in% main_types) %>% 
  ggplot(aes(x = phylum, 
             y = reorder(type, type_norm_abso_proportion_sum), 
             fill = log10(type_norm_abso_proportion_sum))) + 
  geom_raster() + 
  scale_fill_gradient2(low = "#31a354", high = "#e9a3c9", mid = "#ffffbf", midpoint = 2) + 
  coord_fixed() + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = 'bold'),
        axis.text.y = element_text(face = 'bold'),
        legend.title = element_text(angle = 90, size = 8)) + 
  labs(y = "Defense systems (average abundance >= 1%)", x = "Phylum", fill = "log10(average RPKM across all samples)")
#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures_new/26-4-8.type_proportions_in_main_phyla
type_phylum %>% filter(phylum %in% main_phylum & type %in% main_types) %>% 
  ggplot(aes(x = phylum, 
             y = reorder(type, type_norm_abso_proportion_sum), 
             fill = type_norm_abso_proportion_sum/phylum_norm_abso_proportion_sum)) + 
  geom_raster() + 
  scale_fill_gradient2(low = "#31a354", high = "#e9a3c9", mid = "#ffffbf", midpoint = 0.15) + 
  coord_fixed() + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = 'bold'),
        axis.text.y = element_text(face = 'bold'),
        legend.title = element_text(angle = 90, size = 8, vjust = 1)) + 
  labs(y = "Defense systems (average abundance >= 1%)", x = "Phylum", fill = "System relative abundance in the corresponding phylum")
#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures_new/26-4-9.type_abundances_within_main_phyla
# ==============================================================================

# ==============================================================================
# summarize RM genes
rm_genes <- defense_genes %>% filter(type == 'RM')
rm_systems <- rm_genes %>% group_by(sys_id.x, genes_count) %>% summarize(system_gene_count = n(), system = 1)
rm_systems %>% ggplot(aes(x = genes_count)) + 
  geom_histogram() # show the size distribution of RM systems
rm_genes_genus <- rm_genes %>% group_by(subject, time_point, genus) %>% 
  summarize(absolute_proportion_sum = sum(absolute_proportion_of_gene), 
            relative_proportion_sum = sum(relative_proportion_of_gene), 
            normalized_absolute_proportion_sum = sum(normalized_absolute_proportion_of_gene), 
            normalized_relative_proportion_sum = sum(normalized_relative_proportion_of_gene))
rm_genes_genus <- rm_genes_genus[, c(1, 2, 3, 6)] # use normalized relative proportion
rm_genes_genus_matrix <- rm_genes_genus %>% spread(key = genus, value = normalized_absolute_proportion_sum, fill = 0)
rm_genes_genus <- rm_genes_genus_matrix %>% gather(key = 'genus', value = 'normalized_absolute_proportion_sum', 
                                                   colnames(rm_genes_genus_matrix)[3:124])
rm_genes_genus_average <- rm_genes_genus %>% group_by(genus) %>% 
  summarize(average_norm_absolute_sum = sum(normalized_absolute_proportion_sum)/100)
rm_genes_genus <- left_join(rm_genes_genus, rm_genes_genus_average, by = 'genus')
rm_genes_genus$rpkm <- rm_genes_genus$normalized_absolute_proportion_sum * 1e9
rm_genes_genus$average_rpkm <- rm_genes_genus$average_norm_absolute_sum * 1e9
rm_genes_genus %>% 
  filter(genus != 'unclassified') %>%
  filter(average_rpkm >= 5) %>%
  ggplot(aes(x = reorder(genus, -average_rpkm), y = log10(rpkm + 1))) + 
  geom_violin() + 
  geom_point(aes(fill = as.character.numeric_version(subject)), shape = 21, alpha = 0.6, position = position_jitter(width = 0.15)) + 
  geom_point(aes(x = reorder(genus, -average_rpkm), y = log10(average_rpkm+1)), 
             shape = 25, color = 'red', fill = 'red', size = 2) +
  # stat_boxplot(geom = 'errorbar', width = 0.3) +
  # stat_summary(fun = 'mean', geom = 'point', shape = 25, color = 'red', fill = 'red', size = 2) + 
  scale_fill_manual(values = color_31) + 
  theme(legend.position = 'None') + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = 'bold.italic')) + 
  labs(x = 'genus (average RPKM across all samples >= 5)')
#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures_new/26-4-10.RM_genes_abundances_in_genus
# ==============================================================================

# ==============================================================================
# summarize Cas genes' taxonomy
cas_genes <- defense_genes %>% filter(type == 'CasFinder')
cas_systems <- cas_genes %>% group_by(sys_id.x, genes_count) %>% summarize(system_gene_count = n(), system = 1)
cas_systems %>% ggplot(aes(x = genes_count)) + 
  geom_histogram() # show the size distribution of cas systems
cas_genes_genus <- cas_genes %>% group_by(subject, time_point, genus) %>% 
  summarize(absolute_proportion_sum = sum(absolute_proportion_of_gene), 
            relative_proportion_sum = sum(relative_proportion_of_gene), 
            normalized_absolute_proportion_sum = sum(normalized_absolute_proportion_of_gene), 
            normalized_relative_proportion_sum = sum(normalized_relative_proportion_of_gene))
cas_genes_genus <- cas_genes_genus[, c(1, 2, 3, 6)] # use normalized relative proportion
cas_genes_genus_matrix <- cas_genes_genus %>% spread(key = genus, value = normalized_absolute_proportion_sum, fill = 0)
cas_genes_genus <- cas_genes_genus_matrix %>% gather(key = 'genus', value = 'normalized_absolute_proportion_sum', 
                                                     colnames(cas_genes_genus_matrix)[3:108])
cas_genes_genus_average <- cas_genes_genus %>% group_by(genus) %>% 
  summarize(average_nocas_absolute_sum = sum(normalized_absolute_proportion_sum)/100)
cas_genes_genus <- left_join(cas_genes_genus, cas_genes_genus_average, by = 'genus')
cas_genes_genus$rpkm <- cas_genes_genus$normalized_absolute_proportion_sum * 1e9
cas_genes_genus$average_rpkm <- cas_genes_genus$average_nocas_absolute_sum * 1e9
cas_genes_genus %>% 
  filter(genus != 'unclassified') %>%
  filter(average_rpkm >= 5) %>%
  ggplot(aes(x = reorder(genus, -average_rpkm), y = log10(rpkm + 1))) + 
  # geom_violin() + 
  geom_boxplot(width = 0.8, outlier.alpha = 0) +
  geom_point(aes(fill = as.character.numeric_version(subject)), shape = 21, alpha = 0.6, position = position_jitter(width = 0.15)) + 
  geom_point(aes(x = reorder(genus, -average_rpkm), y = log10(average_rpkm+1)), 
             shape = 25, color = 'red', fill = 'red', size = 2) +
  stat_boxplot(geom = 'errorbar', width = 0.3) +
  # stat_summary(fun = 'mean', geom = 'point', shape = 25, color = 'red', fill = 'red', size = 2) + 
  scale_fill_manual(values = color_31) + 
  theme(legend.position = 'None') + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = 'bold.italic')) + 
  labs(x = 'genus (average RPKM across all samples >= 5)')
#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures_new/26-4-11.cas_genes_abundances_in_genus
# ==============================================================================

# ==============================================================================
# summarize paris genes' taxonomy
paris_genes <- defense_genes %>% filter(type == 'Paris')
paris_systems <- paris_genes %>% group_by(sys_id.x, genes_count) %>% summarize(system_gene_count = n(), system = 1)
paris_systems %>% ggplot(aes(x = genes_count)) + 
  geom_histogram() # show the size distribution of paris systems
paris_genes_genus <- paris_genes %>% group_by(subject, time_point, genus) %>% 
  summarize(absolute_proportion_sum = sum(absolute_proportion_of_gene), 
            relative_proportion_sum = sum(relative_proportion_of_gene), 
            normalized_absolute_proportion_sum = sum(normalized_absolute_proportion_of_gene), 
            normalized_relative_proportion_sum = sum(normalized_relative_proportion_of_gene))
paris_genes_genus <- paris_genes_genus[, c(1, 2, 3, 6)] # use normalized relative proportion
paris_genes_genus_matrix <- paris_genes_genus %>% spread(key = genus, value = normalized_absolute_proportion_sum, fill = 0)
paris_genes_genus <- paris_genes_genus_matrix %>% gather(key = 'genus', value = 'normalized_absolute_proportion_sum', 
                                                         colnames(paris_genes_genus_matrix)[3:34])
paris_genes_genus_average <- paris_genes_genus %>% group_by(genus) %>% 
  summarize(average_noparis_absolute_sum = sum(normalized_absolute_proportion_sum)/100)
paris_genes_genus <- left_join(paris_genes_genus, paris_genes_genus_average, by = 'genus')
paris_genes_genus$rpkm <- paris_genes_genus$normalized_absolute_proportion_sum * 1e9
paris_genes_genus$average_rpkm <- paris_genes_genus$average_noparis_absolute_sum * 1e9
paris_genes_genus %>% 
  # filter(genus != 'unclassified') %>%
  filter(average_rpkm >= 0.5) %>%
  ggplot(aes(x = reorder(genus, -average_rpkm), y = log10(rpkm + 1))) + 
  # geom_violin() + 
  geom_boxplot(width = 0.8, outlier.alpha = 0) +
  geom_point(aes(fill = as.character.numeric_version(subject)), shape = 21, alpha = 0.6, position = position_jitter(width = 0.15)) + 
  geom_point(aes(x = reorder(genus, -average_rpkm), y = log10(average_rpkm+1)), 
             shape = 25, color = 'red', fill = 'red', size = 2) +
  stat_boxplot(geom = 'errorbar', width = 0.3) +
  # stat_summary(fun = 'mean', geom = 'point', shape = 25, color = 'red', fill = 'red', size = 2) + 
  scale_fill_manual(values = color_31) + 
  theme(legend.position = 'None') + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = 'bold.italic')) + 
  labs(x = 'genus (average RPKM across all samples >= 0.5)')
#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures_new/26-4-12.paris_genes_abundances_in_genus
# ==============================================================================

# ==============================================================================
# summarize wadjet genes' taxonomy
wadjet_genes <- defense_genes %>% filter(type == 'Wadjet')
wadjet_systems <- wadjet_genes %>% group_by(sys_id.x, genes_count) %>% summarize(system_gene_count = n(), system = 1)
wadjet_systems %>% ggplot(aes(x = genes_count)) + 
  geom_histogram() # show the size distribution of wadjet systems
wadjet_genes_genus <- wadjet_genes %>% group_by(subject, time_point, genus) %>% 
  summarize(absolute_proportion_sum = sum(absolute_proportion_of_gene), 
            relative_proportion_sum = sum(relative_proportion_of_gene), 
            normalized_absolute_proportion_sum = sum(normalized_absolute_proportion_of_gene), 
            normalized_relative_proportion_sum = sum(normalized_relative_proportion_of_gene))
wadjet_genes_genus <- wadjet_genes_genus[, c(1, 2, 3, 6)] # use normalized relative proportion
wadjet_genes_genus_matrix <- wadjet_genes_genus %>% spread(key = genus, value = normalized_absolute_proportion_sum, fill = 0)
wadjet_genes_genus <- wadjet_genes_genus_matrix %>% gather(key = 'genus', value = 'normalized_absolute_proportion_sum', 
                                                           colnames(wadjet_genes_genus_matrix)[3:25])
wadjet_genes_genus_average <- wadjet_genes_genus %>% group_by(genus) %>% 
  summarize(average_nowadjet_absolute_sum = sum(normalized_absolute_proportion_sum)/100)
wadjet_genes_genus <- left_join(wadjet_genes_genus, wadjet_genes_genus_average, by = 'genus')
wadjet_genes_genus$rpkm <- wadjet_genes_genus$normalized_absolute_proportion_sum * 1e9
wadjet_genes_genus$average_rpkm <- wadjet_genes_genus$average_nowadjet_absolute_sum * 1e9
wadjet_genes_genus %>% 
  # filter(genus != 'unclassified') %>%
  filter(average_rpkm >= 0.5) %>%
  ggplot(aes(x = reorder(genus, -average_rpkm), y = log10(rpkm + 1))) + 
  # geom_violin() + 
  geom_boxplot(width = 0.8, outlier.alpha = 0) +
  geom_point(aes(fill = as.character.numeric_version(subject)), shape = 21, alpha = 0.6, position = position_jitter(width = 0.15)) + 
  geom_point(aes(x = reorder(genus, -average_rpkm), y = log10(average_rpkm+1)), 
             shape = 25, color = 'red', fill = 'red', size = 2) +
  stat_boxplot(geom = 'errorbar', width = 0.3) +
  # stat_summary(fun = 'mean', geom = 'point', shape = 25, color = 'red', fill = 'red', size = 2) + 
  scale_fill_manual(values = color_31) + 
  theme(legend.position = 'None') + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = 'bold.italic')) + 
  labs(x = 'genus (average RPKM across all samples >= 0.5)')
#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures_new/26-4-13.wadjet_genes_abundances_in_genus
# ==============================================================================

# ==============================================================================
# summarize abid genes' taxonomy
abid_genes <- defense_genes %>% filter(type == 'AbiD')
abid_systems <- abid_genes %>% group_by(sys_id.x, genes_count) %>% summarize(system_gene_count = n(), system = 1)
abid_systems %>% ggplot(aes(x = genes_count)) + 
  geom_histogram() # show the size distribution of abid systems
abid_genes_genus <- abid_genes %>% group_by(subject, time_point, genus) %>% 
  summarize(absolute_proportion_sum = sum(absolute_proportion_of_gene), 
            relative_proportion_sum = sum(relative_proportion_of_gene), 
            normalized_absolute_proportion_sum = sum(normalized_absolute_proportion_of_gene), 
            normalized_relative_proportion_sum = sum(normalized_relative_proportion_of_gene))
abid_genes_genus <- abid_genes_genus[, c(1, 2, 3, 6)] # use normalized relative proportion
abid_genes_genus_matrix <- abid_genes_genus %>% spread(key = genus, value = normalized_absolute_proportion_sum, fill = 0)
abid_genes_genus <- abid_genes_genus_matrix %>% gather(key = 'genus', value = 'normalized_absolute_proportion_sum', 
                                                       colnames(abid_genes_genus_matrix)[3:72])
abid_genes_genus_average <- abid_genes_genus %>% group_by(genus) %>% 
  summarize(average_noabid_absolute_sum = sum(normalized_absolute_proportion_sum)/100)
abid_genes_genus <- left_join(abid_genes_genus, abid_genes_genus_average, by = 'genus')
abid_genes_genus$rpkm <- abid_genes_genus$normalized_absolute_proportion_sum * 1e9
abid_genes_genus$average_rpkm <- abid_genes_genus$average_noabid_absolute_sum * 1e9
abid_genes_genus %>% 
  # filter(genus != 'unclassified') %>%
  filter(average_rpkm >= 1) %>%
  ggplot(aes(x = reorder(genus, -average_rpkm), y = log10(rpkm + 1))) + 
  # geom_violin() + 
  geom_boxplot(width = 0.8, outlier.alpha = 0) +
  geom_point(aes(fill = as.character.numeric_version(subject)), shape = 21, alpha = 0.6, position = position_jitter(width = 0.15)) + 
  geom_point(aes(x = reorder(genus, -average_rpkm), y = log10(average_rpkm+1)), 
             shape = 25, color = 'red', fill = 'red', size = 2) +
  stat_boxplot(geom = 'errorbar', width = 0.3) +
  # stat_summary(fun = 'mean', geom = 'point', shape = 25, color = 'red', fill = 'red', size = 2) + 
  scale_fill_manual(values = color_31) + 
  theme(legend.position = 'None') + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = 'bold.italic')) + 
  labs(x = 'genus (average RPKM across all samples >= 1)')
#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures_new/26-4-14.abid_genes_abundances_in_genus
# ==============================================================================

# ==============================================================================
# summarize the taxonomic distributions (genus) of different types of defense systems, with normalized absolute proportion of the genes (RPKM)
type_genus <- defense_genes %>% group_by(type, genus) %>% 
  summarize(type_norm_abso_proportion_sum = sum(normalized_absolute_proportion_of_gene))
type_genus <- drop_na(type_genus)
type_genus$rpkm <- type_genus$type_norm_abso_proportion_sum * 1e9
type_genus <- type_genus[, -3]
type_genus_matrix <- type_genus %>% spread(key = genus, value = rpkm, fill = 0)
type_genus <- type_genus_matrix %>% gather(key = 'genus', value = 'rpkm', colnames(type_genus_matrix)[2:187])
genus_rpkm_sum <- type_genus %>% group_by(genus) %>% summarize(genus_rpkm_sum = sum(rpkm))
type_rpkm_sum <- type_genus %>% group_by(type) %>% summarize(type_rpkm_sum = sum(rpkm))
type_genus <- left_join(type_genus, genus_rpkm_sum, by = 'genus')
type_genus <- left_join(type_genus, type_rpkm_sum, by = 'type')

# use previously summarized taxonomic composition table to summarize the top-ranked genus
total_genus <- read.delim('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/20.taxonomy_profiling/mmseqs2+bowtie2/mmseqs2+bowtie2_summarized_reports/20-5-3.bowtie2_genus_report.txt', 
                          sep = '\t', header = T)
genus_total_counts <- total_genus %>% group_by(genus) %>% summarize(genus_total_counts = sum(genus_count))
main_genus <- genus_total_counts[genus_total_counts$genus_total_counts >= 8e6, ][['genus']]
main_types <- unique(defense_types_norm_relative[defense_types_norm_relative$average_type_abundance >= 0.01,][['type']])

type_genus <- left_join(type_genus, genus_total_counts, by = 'genus')
type_genus %>% 
  filter(type %in% main_types) %>%
  filter(genus %in% main_genus & genus != 'unclassified') %>% 
  ggplot(aes(x = reorder(genus, -genus_total_counts), 
             y = reorder(type, type_rpkm_sum), 
             fill = log10(rpkm))) + 
  geom_raster() + 
  scale_fill_gradient2(low = "#31a354", high = "#e9a3c9", mid = "#ffffbf", midpoint = 1.5) + 
  coord_fixed() + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = 'bold.italic'),
        axis.text.y = element_text(face = 'bold'),
        legend.title = element_text(angle = 90, size = 8)) + 
  labs(y = "Defense systems (average abundance >= 1%)", x = "Genera of top abundances", fill = "log10(average RPKM across all samples)")
#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures_new/26-4-15.type_proportions_in_main_genera
type_genus %>% 
  filter(type %in% main_types) %>% 
  filter(genus %in% main_genus & genus != 'unclassified') %>% 
  ggplot(aes(x = reorder(genus, -genus_total_counts), 
             y = reorder(type, type_rpkm_sum),
             fill = rpkm/genus_rpkm_sum)) + 
  geom_raster() + 
  scale_fill_gradient2(low = "#31a354", high = "#e9a3c9", mid = "#ffffbf", midpoint = 0.25) + 
  coord_fixed() + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = 'bold.italic'),
        axis.text.y = element_text(face = 'bold'),
        legend.title = element_text(angle = 90, size = 8)) + 
  labs(y = "Defense systems (average abundance >= 1%)", x = "Genera of top abundances", fill = "System relative abundance in the corresponding genus")
#~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures_new/26-4-16.type_abundances_within_main_genera
# ==============================================================================