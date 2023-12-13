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

#

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