# abundance analysis of defense systems are not normalized in this script.

library(randomcoloR)
library(tidyr)
library(dplyr)
library(stringr)

getwd()
setwd("/home/zhanggaopu/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts")

# summarize genes/systems predicted with DefenseFinder, qualitative, not quantitative.
genes <- read.table('../intermediates/25.defensefinder/defensefinder/merged_defensefinder_genes.tsv', 
                    header = T, sep = '\t')
systems <- read.table('../intermediates/25.defensefinder/defensefinder/merged_defensefinder_systems.tsv', 
                      header = T, sep = '\t')

genes <-cbind(genes$hit_id, genes %>% separate(hit_id, sep = '_', into = c('individual', 'k141', 'contig_code', 'flag', 'multi', 'contig_len', 'CDS_code')))
systems <- cbind(systems$sys_id, systems %>% separate(sys_id, sep = '_', into = c('individual')))

systems_summary <- systems %>% group_by(individual, type) %>% 
  summarise(system_count = n())
n_distinct(systems_summary$type)
systems_summary %>% ggplot(aes(x = individual, y = system_count)) + 
  geom_col(aes(fill = type))
individual_systems_count <- systems %>% group_by(individual) %>% 
  summarise(individual_system_count = n())
systems_summary <- left_join(systems_summary, individual_systems_count, by = 'individual')
systems_summary <- systems_summary %>% mutate(system_proportion = system_count/individual_system_count)
systems_summary <- systems_summary %>% mutate(system_legend = case_when(
  system_proportion >= 0.01 ~ type, 
  system_proportion < 0.01 ~ ' others'
))
n_distinct(systems_summary$system_legend)
color_28 <- distinctColorPalette(27)
color_28 <- c('#bdbdbd', color_28)
systems_summary %>% ggplot(aes(x = individual, y = system_proportion)) + 
  geom_col(aes(fill = system_legend)) + 
  scale_fill_manual(values = color_28)

systems_summary %>% filter(system_legend != 'CasFinder' & system_legend != 'RM') %>%
  ggplot(aes(x = individual, y = system_proportion)) + 
  geom_col(aes(fill = system_legend)) + 
  scale_fill_manual(values = color_28) + 
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))

systems_summary %>% ggplot(aes(x = reorder(system_legend, -system_proportion), y = system_proportion)) + 
  geom_col(aes(fill = type))

#-------------------------------------------------------------------------------
# Summarize the anti-phage genes and their corresponding systems.
systems <- read.delim('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/25.defensefinder/defensefinder/merged_defensefinder_systems.tsv', 
                           sep = '\t', header = T)
total_orfs <- str_split_fixed(systems$protein_in_syst, ',', 16)
systems <- cbind(systems, total_orfs)
systems[systems == ''] <- NA
systems %>% ggplot(aes(x = genes_count)) +
  geom_histogram(binwidth = 0.5)


cas_systems <- systems %>% filter(type == "CasFinder")
cas_genes <- str_split_fixed(cas_systems$name_of_profiles_in_sys, ',', 16) # split the genes of each system into separated columns
cas_genes[cas_genes == ''] <- NA
cas_systems <- cbind(cas_systems, cas_genes)
cas_systems %>% ggplot(aes(x = genes_count)) + 
  geom_histogram(binwidth = 0.5)
cas_systems <- cas_systems %>% gather(key = 'cas_gene_count', value = 'cas_gene', `1`, `2`, `3`, `4`, 
                                      `5`, `6`, `7`, `8`, `9`, `10`, `11`, `12`, `13`, `14`, `15`, `16`)

#-------------------------------------------------------------------------------
# Summarize the counts of genes from featureCounts.
counts <- read.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/featureCounts/merged_gene_count.txt', 
                     sep = ' ', header = F)
colnames(counts) <- c('gene_CDS', '10', '1', '2', '3', '4', '5', '6', '7',
                      '8', '9', 'individual')
counts$gene_CDS <- gsub(pattern = '\\^', replacement = '=', x = counts$gene_CDS) # ^ is special character.
counts <- counts %>% separate(gene_CDS, sep = '[+]', into = c('gene', 'CDS')) # sep="" requires regular expression!
counts %>% write.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/featureCounts/26-4.summarized_merged_gene_count.txt', 
                       sep = '\t', quote = F, row.names = F, col.names = T)

antiphage_counts <- counts %>% filter(gene != "non-antiphage-gene") 
antiphage_counts <- antiphage_counts %>% gather(key = 'time_point', value = 'mapping_count', 
                                                `1`, `2`, `3`, `4`, `5`, 
                                                `6`, `7`, `8`, `9`, `10`)

# annotation_list <- read.delim('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/26-1.total_CDSs_annotation_list.txt', 
                              # sep = '\t', header = T)
# annotation_list$CDS <- gsub(pattern = "\\^", replacement = "=", x = annotation_list$CDS)
# antiphage_counts <- left_join(antiphage_counts, annotation_list, by = "CDS")

total_antiphage_systems <- read.delim('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/25.defensefinder/defensefinder/merged_defensefinder_systems.tsv', 
                                      sep = '\t', header = T)
total_antiphage_CDS <- str_split_fixed(total_antiphage_systems$protein_in_syst, ',', 16)
total_antiphage_systems <- cbind(total_antiphage_systems, total_antiphage_CDS)
total_antiphage_systems[total_antiphage_systems == ''] <- NA
total_antiphage_genes <- total_antiphage_systems %>% gather(key = 'gene_code_in_sys', value = 'CDS', colnames(total_antiphage_systems)[9:24])
total_antiphage_genes <- total_antiphage_genes[!(is.na(total_antiphage_genes$CDS)), ]

n_distinct(total_antiphage_genes$CDS)
duplicated_genes <- total_antiphage_genes[duplicated(total_antiphage_genes$CDS), c(1, 7, 10)]
# some predicted genes are shared by two or more systems. The counts of these genes were repeatedly counted. May need further de-complexing.

antiphage_counts <- left_join(antiphage_counts, total_antiphage_genes, by = 'CDS')
antiphage_counts$individual <- as.numeric(antiphage_counts$individual)
antiphage_counts$time_point <- as.numeric(antiphage_counts$time_point)
individual_total_counts <- antiphage_counts %>% group_by(individual, time_point) %>% 
  summarize(total_counts = sum(mapping_count))

type_level <- antiphage_counts %>% group_by(individual, time_point, type) %>% 
  summarize(count = sum(mapping_count))
type_level <- left_join(type_level, individual_total_counts, by = c('individual', 'time_point'))
type_level <- type_level %>% mutate(type_proportion = count/total_counts)
type_level <- type_level %>% mutate(legend = case_when(
  type_proportion >= 0.02 ~ type,
  type_proportion < 0.02 ~ 'others'
))
type_level$legend <- reorder(type_level$legend, type_level$type_proportion)
n_distinct(type_level$legend)
color_38 <- distinctColorPalette(37)
color_43 <- c("#F135CC", "#E1E4CF", "#A7AFDF", "#E596E2", "#86EB2D", "#DBE5A5", 
              "#DFA4A0", "#95D9E6", "#A681DE", "#E7509A", "#5563DA", "#B9EAC9", 
              "#67E6C4", "#C976E8", "#5DE55F", "#E666DC", "#D7D569", "#D36C51", 
              "#8E9DEC", "#D5E3F0", "#63E79C", "#9A3C8F", "#9DE674", "#68A594",
              "#565DA2", "#DAC499", "#ACE599", "#E43E53", "#5EA9DA", "#DCC5D9",
              "#CD2FED", "#894BDB", "#C89A61", "#69E7E7", "#E7AFDA", "#DDEA4E", 
              "#6B9B4F", "#837F84", "#633BEE", "#D47497", "#E6A637", "#bdbdbd", 
              "#737373")
color_38 <- c('#bdbdbd', color_43[1:37])
type_level %>% ggplot(aes(x = time_point, y = type_proportion)) + 
  geom_col(aes(fill = legend)) + 
  facet_wrap(individual ~ ., nrow = 2) + 
  scale_fill_manual(values = color_38) +
  scale_x_continuous(limits = c(0.5, 10.5), breaks = c(1, 2, 3, 4, 5,
                                                       6, 7, 8, 9, 10)) + 
  labs(fill = 'Defense system type')
# ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures/26-4-1.defense_system_type_total_summary

type_level %>% filter(legend != 'RM' & legend != 'CasFinder') %>% ggplot(aes(x = time_point, y = type_proportion)) + 
  geom_col(aes(fill = legend)) + 
  facet_wrap(individual ~ ., nrow = 2) + 
  scale_fill_manual(values = color_38) +
  scale_x_continuous(limits = c(0.5, 10.5), breaks = c(1, 2, 3, 4, 5,
                                                       6, 7, 8, 9, 10)) + 
  labs(fill = 'Defense system type')
# ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures/26-4-2.defense_system_type_no_Cas_RM_summary

type_average <- type_level %>% group_by(type) %>% 
  summarise(type_average_proportion = sum(type_proportion)/100)

type_level <- left_join(type_level, type_average, by = 'type')

type_level %>% 
  filter(type_average_proportion >= 0.01) %>%
  ggplot(aes(x = reorder(type, -type_average_proportion), y = type_proportion)) + 
  geom_boxplot(outlier.alpha = 0, notch = T, notchwidth = 0.8) + 
  geom_point(aes(fill = as.character.numeric_version(individual)), shape = 21, 
             alpha = 0.7, position = position_jitter(0.1)) + 
  scale_fill_manual(values = color_20) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
# ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures/26-4-5.defense_system_type_average_proportion

#-------------------------------------------------------------------------------
# summarize Cas genes
cas_counts <- antiphage_counts %>% filter(type == 'CasFinder')
cas_genes <- separate(cas_counts, col = 'gene', sep = '_', into = c("cas_gene"))
cas_individual_total_counts <- cas_genes %>% group_by(individual, time_point) %>% 
  summarize(cas_total_counts = sum(mapping_count))
cas_gene_level <- cas_genes %>% group_by(individual, time_point, cas_gene) %>% 
  summarize(cas_gene_count = sum(mapping_count))
cas_gene_level <- left_join(cas_gene_level, cas_individual_total_counts, by = c('individual', 'time_point'))
cas_gene_level <- cas_gene_level %>% mutate(cas_gene_proportion = cas_gene_count/cas_total_counts)
cas_gene_level <- cas_gene_level %>% mutate(legend = case_when(
  cas_gene_proportion >= 0.05 ~ cas_gene,
  cas_gene_proportion < 0.05 ~ 'others'
))

cas_gene_average <- cas_gene_level %>% group_by(cas_gene) %>% 
  summarise(average_gene_proportion = sum(cas_gene_proportion)/100)
cas_gene_level <- left_join(cas_gene_level, cas_gene_average, by = 'cas_gene')
cas_gene_level <- cas_gene_level %>% mutate(legend_average = case_when(
  average_gene_proportion >= 0.04 ~ cas_gene,
  average_gene_proportion < 0.04 ~ 'others'
))

cas_gene_level$legend_average <- reorder(cas_gene_level$legend_average, cas_gene_level$average_gene_proportion)
n_distinct(cas_gene_level$legend_average)
color_17 <- c('#bdbdbd', color_43[1:16])
cas_gene_level$legend <- reorder(cas_gene_level$legend, cas_gene_level$cas_gene_proportion)
cas_gene_level %>% ggplot(aes(x = time_point, y = cas_gene_proportion)) + 
  geom_col(aes(fill = legend_average)) + 
  scale_fill_manual(values = color_34) + 
  facet_wrap(individual~., nrow = 2)

cas_gene_level$legend_average <- reorder(cas_gene_level$legend_average, -cas_gene_level$average_gene_proportion)
cas_gene_level %>% filter(legend_average != 'others') %>%
  ggplot(aes(x = legend_average, y = cas_gene_proportion)) + 
  geom_boxplot(width = 0.6, outlier.colour = NA, notch = TRUE, notchwidth = 0.8) + 
  geom_point(aes(fill = legend_average), 
             shape = 21, size = 2,
             position = position_jitter(width = 0.1)) + 
  scale_fill_manual(values = color_43) + 
  theme(axis.text.x = element_text(angle = 10, size = 12, vjust = 0.8)) + 
  guides(fill = "none")
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# summarize the taxonomic source of cas genes.
total_taxon <- read.delim('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/15.mmseqs_taxonomy/taxonomy_new/mmseqs2_total_summary.txt', 
                          sep = '\t', header = T)
antiphage_contigs <- as.data.frame(antiphage_counts$CDS)
colnames(antiphage_contigs)[1] <- 'CDS'
antiphage_contigs <- separate(antiphage_contigs, col = CDS, sep = "len=", into = c('V1', 'V2'))
antiphage_contigs <- separate(antiphage_contigs, col = V2, sep = "_", into = c('V2', 'V3'))
antiphage_contigs$contig <- paste(antiphage_contigs$V1, antiphage_contigs$V2, sep = "len=")

antiphage_counts$contig <- antiphage_contigs$contig
antiphage_counts <- left_join(antiphage_counts, total_taxon, by = 'contig')

type_level_genus <- antiphage_counts %>% group_by(individual, time_point, type, genus) %>% 
  summarize(type_genus_count = sum(mapping_count))
type_level_genus <- left_join(type_level_genus, individual_total_counts, by = c('individual', 'time_point')) 
type_level_genus <- left_join(type_level_genus, type_level[, c(1, 2, 3, 4)], by = c('individual', 'time_point', 'type'))
type_level_genus <- type_level_genus %>% mutate(genus_proportion_in_type = type_genus_count/count)
cas_genus <- type_level_genus %>% filter(type == 'CasFinder')
cas_genus <- cas_genus %>% mutate(genus_legend = case_when(
  genus_proportion_in_type >= 0.03 ~ genus,
  genus_proportion_in_type < 0.03 ~ 'others'
))
n_distinct(cas_genus$genus_legend)
#cas_genus$genus_legend <- reorder(cas_genus$genus_legend, cas_genus$genus_proportion_in_type)
cas_genus %>% ggplot(aes(x = time_point, y = genus_proportion_in_type)) + 
  geom_col(aes(fill = genus_legend)) + 
  scale_fill_manual(values = color_43) + 
  facet_wrap(individual~., nrow = 1)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# compare taxonomic feature of cas and crispr
# read the table of CRISPR taxonomic composition
total_reads_stat <- read.delim('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/19.bowtie2_reads_contig/figures/19-3-1.crispr_reads_taxonomy_stat.txt',
                               header = T, sep = '\t')
sample_read_sum <- total_reads_stat %>% group_by(sample_info) %>% 
  summarise(read_count = n())
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
crispr_taxon <- total_reads_genus_stat[, c(1, 2, 4, 5)]

cas_taxon <- cas_genus[, c(1, 2, 4, 8)]
colnames(crispr_taxon)[4] <- 'genus_proportion_CRISPR'
colnames(cas_taxon)[4] <- 'genus_proportion_Cas'
crispr_cas_taxon <- full_join(cas_taxon, crispr_taxon, by = c('individual', 'time_point', 'genus'))
crispr_cas_taxon <- crispr_cas_taxon[str_detect(crispr_cas_taxon$genus, 'g__'), ]
crispr_cas_taxon <- crispr_cas_taxon[!(is.na(crispr_cas_taxon$genus)), ]
crispr_cas_taxon[is.na(crispr_cas_taxon)] <- 0
crispr_cas_taxon %>% 
  filter(genus_proportion_Cas >= 0.005 & genus_proportion_CRISPR >= 0.005) %>%
  ggplot(aes(x = genus_proportion_CRISPR, y = genus_proportion_Cas)) + 
  geom_point() + 
  geom_smooth(method = 'lm') + 
  scale_x_continuous(limits = c(0, 0.45)) + 
  scale_y_continuous(limits = c(0, 0.45)) + 
  geom_vline(xintercept = 0.005, linetype = 'dashed', color = '#bdbdbd') + 
  geom_hline(yintercept = 0.005, linetype = 'dashed', color = '#bdbdbd') + 
  geom_abline(slope = 1, linetype = 'dashed', color = '#bdbdbd')
cas_crispr_lm <- crispr_cas_taxon %>% filter(genus_proportion_Cas >= 0.005 & genus_proportion_CRISPR >= 0.005) %>% 
  lm(formula = genus_proportion_Cas ~ genus_proportion_CRISPR)
summary(cas_crispr_lm)
# "~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures/26-4-3.crispr_cas_taxon_association"

crispr_cas_taxon %>%
  # filter(genus_proportion_Cas >= 0.005 & genus_proportion_CRISPR >= 0.005) %>%
  filter(genus %in% c('g__Bacteroides', 'g__Escherichia', 'g__Ligilactobacillus', 'g__Megamonas', 'g__Ruminococcus')) %>% 
  ggplot(aes(x = genus_proportion_CRISPR, y = genus_proportion_Cas)) + 
  geom_point(aes(fill = genus), shape = 21, size = 2) + 
  geom_smooth(method = 'lm') +
  scale_x_continuous(limits = c(0, NA)) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_manual(values = color_20[c(2, 9, 12, 13, 19)]) + 
  facet_wrap(.~ genus, scales = "free", nrow = 2) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) + 
  theme(legend.position = "None") + 
  theme(strip.text = element_text(face = 'italic'))
# "~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures/26-4-4.crispr_cas_taxon_association_main_genus"

crispr_cas_taxon_tidy <- crispr_cas_taxon %>% gather(key = 'object', value = 'proportion', 
                                                     genus_proportion_Cas, genus_proportion_CRISPR)
crispr_cas_taxon_tidy$sample_info <- paste(crispr_cas_taxon_tidy$individual, crispr_cas_taxon_tidy$time_point, sep = '_')
crispr_cas_taxon_tidy %>% 
  filter(genus %in% c('g__Bacteroides', 'g__Escherichia', 'g__Ligilactobacillus', 'g__Megamonas', 'g__Ruminococcus')) %>% 
  ggplot(aes(x = object, y = proportion)) + 
  geom_point(aes(fill = genus), shape = 21) + 
  geom_line(aes(group = sample_info, color = genus)) + 
  facet_grid(. ~ genus)

# ------------------------------------------------------------------------------
# summarize taxonomic soure of Wadjet.
wadjet_counts <- antiphage_counts %>% filter(type == 'Wadjet')
wadjet_individual_sum <- wadjet_counts %>% group_by(individual, time_point) %>% 
  summarize(individual_wadjet_counts = sum(mapping_count))
wadjet_counts <- left_join(wadjet_counts, wadjet_individual_sum, by = c('individual', 'time_point'))
wadjet_genus <- wadjet_counts %>% group_by(individual, time_point, type, genus) %>% 
  summarize(genus_count_wadjet = sum(mapping_count))
wadjet_genus <- left_join(wadjet_genus, wadjet_individual_sum, by = c('individual', 'time_point'))
wadjet_genus <- wadjet_genus %>% mutate(genus_proportion_wadjet = genus_count_wadjet/individual_wadjet_counts)
n_distinct(wadjet_genus$genus)
wadjet_genus %>% ggplot(aes(x = time_point, y = genus_proportion_wadjet)) + 
  geom_col(aes(fill = reorder(genus, genus_proportion_wadjet))) + 
  facet_wrap(individual ~ ., nrow = 2) + 
  scale_fill_manual(values = color_34)
# "~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures/26-4-7.wadjet_genus_proportion"


wadjet_genus_average <- wadjet_genus %>% group_by(genus) %>% 
  summarize(average_proportion = sum(genus_proportion_wadjet)/100)
wadjet_genus <- full_join(wadjet_genus, wadjet_genus_average, by = 'genus')
wadjet_genus %>%
  filter(genus != 'unclassified' & average_proportion >= 0.01) %>% 
  ggplot(aes(x = reorder(genus, -average_proportion), y = genus_proportion_wadjet)) + 
  geom_boxplot(outlier.alpha = 0) +
  geom_point(aes(fill = as.character.numeric_version(individual)), 
             shape = 21, position = position_jitter(0.1), alpha = 0.5) +
  scale_fill_manual(values = color_20,) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = 'italic'))
# "~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures/26-4-6.wadjet_genus_proportion"


# summarize taxonomic source of BREX
brex_counts <- antiphage_counts %>% filter(type == 'BREX')
brex_individual_sum <- brex_counts %>% group_by(individual, time_point) %>% 
  summarize(individual_brex_counts = sum(mapping_count))
brex_counts <- left_join(brex_counts, brex_individual_sum, by = c('individual', 'time_point'))
brex_genus <- brex_counts %>% group_by(individual, time_point, type, genus) %>% 
  summarize(genus_count_brex = sum(mapping_count))
brex_genus <- left_join(brex_genus, brex_individual_sum, by = c('individual', 'time_point'))
brex_genus <- brex_genus %>% mutate(genus_proportion_brex = genus_count_brex/individual_brex_counts)
n_distinct(brex_genus$genus)
brex_genus %>% ggplot(aes(x = time_point, y = genus_proportion_brex)) + 
  geom_col(aes(fill = reorder(genus, genus_proportion_brex))) + 
  facet_wrap(individual ~ ., nrow = 2) + 
  scale_fill_manual(values = color_34)
# "~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures/26-4-8.brex_genus_proportion"


brex_genus_average <- brex_genus %>% group_by(genus) %>% 
  summarize(average_proportion = sum(genus_proportion_brex)/100)
brex_genus <- full_join(brex_genus, brex_genus_average, by = 'genus')
brex_genus %>%
  filter(genus != 'unclassified' & average_proportion >= 0.01) %>% 
  ggplot(aes(x = reorder(genus, -average_proportion), y = genus_proportion_brex)) + 
  geom_boxplot(outlier.alpha = 0) +
  geom_point(aes(fill = as.character.numeric_version(individual)), 
             shape = 21, position = position_jitter(0.1), alpha = 0.5) +
  scale_fill_manual(values = color_20,) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = 'italic'))
# "~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures/26-4-9.brex_genus_proportion"


# summarize taxonomic source of gabija
gabija_counts <- antiphage_counts %>% filter(type == 'Gabija')
gabija_individual_sum <- gabija_counts %>% group_by(individual, time_point) %>% 
  summarize(individual_gabija_counts = sum(mapping_count))
gabija_counts <- left_join(gabija_counts, gabija_individual_sum, by = c('individual', 'time_point'))
gabija_genus <- gabija_counts %>% group_by(individual, time_point, type, genus) %>% 
  summarize(genus_count_gabija = sum(mapping_count))
gabija_genus <- left_join(gabija_genus, gabija_individual_sum, by = c('individual', 'time_point'))
gabija_genus <- gabija_genus %>% mutate(genus_proportion_gabija = genus_count_gabija/individual_gabija_counts)
n_distinct(gabija_genus$genus)
gabija_genus %>% ggplot(aes(x = time_point, y = genus_proportion_gabija)) + 
  geom_col(aes(fill = reorder(genus, genus_proportion_gabija))) + 
  facet_wrap(individual ~ ., nrow = 2) + 
  scale_fill_manual(values = color_34)

gabija_genus_average <- gabija_genus %>% group_by(genus) %>% 
  summarize(average_proportion = sum(genus_proportion_gabija)/100)
gabija_genus <- full_join(gabija_genus, gabija_genus_average, by = 'genus')
gabija_genus %>%
  filter(genus != 'unclassified' & average_proportion >= 0.01) %>% 
  ggplot(aes(x = reorder(genus, -average_proportion), y = genus_proportion_gabija)) + 
  geom_boxplot(outlier.alpha = 0) +
  geom_point(aes(fill = as.character.numeric_version(individual)), 
             shape = 21, position = position_jitter(0.1), alpha = 0.5) +
  scale_fill_manual(values = color_20,) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, face = 'italic'))

# summarize g__Bacteroides's machineries
bacteroides_counts <- antiphage_counts %>% filter(genus == 'g__Bacteroides')
bacteroides_types <- bacteroides_counts %>% group_by(individual, time_point, type) %>% 
  summarize(type_counts = sum(mapping_count))
individual_bacteroides_sum <- bacteroides_counts %>% group_by(individual, time_point) %>% 
  summarize(bacteroides_counts = sum(mapping_count))
bacteroides_types <- left_join(bacteroides_types, individual_bacteroides_sum, by = c('individual', 'time_point'))
bacteroides_types <- bacteroides_types %>% mutate(type_proportion = type_counts/bacteroides_counts)
bacteroides_types <- bacteroides_types %>% mutate(legend = case_when(
  type_proportion >= 0.05 ~ type,
  type_proportion < 0.05 ~ 'others'
))
n_distinct(bacteroides_types$legend)
bacteroides_types$legend <- reorder(bacteroides_types$legend, bacteroides_types$type_proportion)
bacteroides_types %>% ggplot(aes(x = time_point, y = type_proportion)) + 
  facet_wrap(individual ~., nrow = 2) +
  geom_col(aes(fill = legend)) + 
  scale_fill_manual(values = color_34)
# "~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures/26-4-10.Bacteroides_machinery_summary"
bacteroides_types_average <- bacteroides_types %>% group_by(type) %>% 
  summarize(type_average_proportion = sum(type_proportion)/100)
bacteroides_types <- left_join(bacteroides_types, bacteroides_types_average, by = 'type')
bacteroides_types %>% 
  filter(type_average_proportion >= 0.02) %>%  
  ggplot(aes(x = reorder(type, -type_average_proportion), y = type_proportion)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_point(aes(fill = as.character.numeric_version(individual)), shape = 21, 
             position = position_jitter(0.2), alpha = 0.7) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
# "~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/figures/26-4-11.Bacteroides_machinery_summary"
