# This script is for summarizing the taxonomic assignments of each contig, which will be used in future analysis of metagenomic and CRISPR taxonomic analysis.
#----------------------------------------------------------------------------------------------------------------------------------------

getwd()

total <- read.csv('../intermediates/15.mmseqs_taxonomy/taxonomy_new/total.csv', 
                    header = F, sep = ';')
# Preprocessing of the data frame.
contig_info <- as.data.frame(total$V1)
contig_info <- contig_info %>% separate(`total$V1`,c("individual", "contig_info_1", "contig_info_2", "flag", "multi", "len"), sep = '_')
individual <- contig_info$individual
contig_len <- contig_info$len
contig_info <- paste(contig_info$contig_info_1, contig_info$contig_info_2, sep = '_')
total$individual <- individual
total$contig_info <- contig_info
total$contig_len <- contig_len
colnames(total)[1:6] <- c("contig", "taxid", "total_orf", "assigned_orf", "corresponding_orf", "-logE")

# Extract taxid
taxid <- total$taxid
write.table(taxid, '../intermediates/15.mmseqs_taxonomy/taxonomy_new/taxid.txt', 
            quote = F, sep = '\t', row.names = F, col.names = F)
# Corresponding taxonomy of the taxid is assigned with TaxonKit.
taxonomy <- read.table('../intermediates/15.mmseqs_taxonomy/taxonomy_new/taxonomy.txt', header = F, sep = '\t', fill = NA)
taxonomy[taxonomy==""] <- 'unclassified'
for (i in 2:ncol(taxonomy)) {
  taxonomy[, i][str_detect(taxonomy[, i], 'unclassified')] <- 'unclassified'
  print(i)
}
taxonomy <- taxonomy[, -2]
colnames(taxonomy) <- c('taxid', 'kingdom', 'phylum', 'class', 'order', 
                        'family', 'genus', 'species')

# Merging the taxonomy information into the total dataframe.
total <- cbind(total, taxonomy[, 2:8])
write.table(total, '../intermediates/15.mmseqs_taxonomy/taxonomy_new/mmseqs2_total_summary.txt', 
          sep = '\t', quote = F, row.names = F, col.names = T)

total <- read.table('../intermediates/15.mmseqs_taxonomy/taxonomy_new/mmseqs2_total_summary.txt',
                    sep = '\t', header = T)

# Taxonomy statistics of each individual
phylum <- total %>% group_by(individual, phylum) %>% 
  summarise(phylum_count = n())
phylum_sum <- phylum %>% group_by(individual) %>% summarise(contig_count = sum(phylum_count))
phylum <- left_join(phylum, phylum_sum)

four_phylum <- c("p__Actinomycetota", "p__Bacillota", "p__Bacteroidota", "p__Pseudomonadota")

"p__sss" != "unclassified" && !("p__sss" %in% four_phylum)

phylum <- phylum %>% mutate(phylum_simplified = case_when(
  phylum == 'unclassified' ~ "unclassified",
  phylum %in% four_phylum ~ phylum,
  !(phylum %in% four_phylum) && phylum != 'unclassified' ~ "others",
))

phylum <- phylum %>% mutate(phylum_type = case_when(
  phylum == 'unclassified' ~ 'unclassified', 
  phylum_count/contig_count >= 0.005 ~ phylum, 
  phylum_count/contig_count < 0.005 ~ 'others'
))

phylum$phylum_simplified <- factor(phylum$phylum_simplified, 
                                   levels = c("unclassified", "others", "p__Pseudomonadota", "p__Actinomycetota", "p__Bacteroidota", "p__Bacillota"))
phylum$individual <- factor(phylum$individual, 
                            levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
phylum %>% ggplot(aes(x = individual, y = phylum_count/contig_count)) + 
  geom_col(aes(fill = phylum_simplified)) + 
  scale_fill_manual(values = c('#969696', '#d9d9d9', '#7fc97f', '#e7298a', '#ffed6f', '#2171b5')) + 
  theme_minimal()

genus <- total %>% group_by(individual, genus) %>% 
  summarise(genus_count = n())
genus_sum <- genus %>% group_by(individual) %>% summarise(contig_count = sum(genus_count))
genus <- left_join(genus, genus_sum)
genus <- genus %>% mutate(genus_simplified = case_when(
  genus == 'unclassified' ~ "unclassified",
  genus_count/contig_count >= 0.01 ~ genus,
  genus_count/contig_count < 0.01 ~ 'others'
))
genus$individual <- factor(genus$individual, 
                            levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
manual_color=c('#8dd3c7','#f16913','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#2171b5','#bc80bd','#ccebc5','#ffed6f',
               '#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#bdbdbd', '#7fc97f','#beaed4','#fdc086', '#d9d9d9', '#969696')
genus %>% ggplot(aes(x = individual, y = genus_count/contig_count)) + 
  geom_col(aes(fill = genus_simplified)) + 
  scale_fill_manual(values = manual_color)+
  theme_minimal()
genus %>% filter(genus_simplified != 'unclassified') %>% ggplot(aes(x = individual, y = genus_count/contig_count)) + 
  geom_col(aes(fill = genus_simplified)) + 
  scale_fill_manual(values = manual_color) + 
  theme_minimal()
