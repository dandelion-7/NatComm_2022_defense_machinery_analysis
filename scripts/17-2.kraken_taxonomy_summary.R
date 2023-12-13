library(stringr)
library(tidyr)
library(dplyr)
getwd()

# Read the total_output.txt file from kraken.
total <- read.csv('../intermediates/17.kraken_taxonomy/kraken_results_summary/kraken_total_output.csv', 
                  header = F, sep = ',')
total <- total %>% separate(V2, c("individual", "time_point", "sample", "group", 
                                  "group_number","repeat", "SRR", "read"))
colnames(total)[c(1, 10, 11, 12, 13)] <- c('classified', 'taxonomy', 'taxid', 
                                           'read_length', 'kmer_mapping')
total <- subset(total, select = -group)
taxid <- total$taxid
write.table(x = taxid, file = '../intermediates/17.kraken_taxonomy/kraken_results_summary/taxid.txt', 
            sep = '\t', quote = F, row.names = F, col.names = F)

# with the taxid.txt file, taxonkit is used to generate the corresponding taxonomical names.
lineage <- read.csv ('../intermediates/17.kraken_taxonomy/kraken_results_summary/lineage.txt', 
                      header = F, sep = '\t')
colnames(lineage) <- c('taxid', 'taxonomy', 'kingdom', 'phylum', 'class', 
                       'order', 'family', 'genus', 'species')
lineage[lineage == ""] <- 'unclassified'
lineage <- lineage[-2]

# Combine the taxonomy information with the CRISPR information
total <- total %>% cbind(lineage)
total <- total[-10]

# Get the number of reads of each CRISPR array (group)
group_sum <- total %>% group_by(individual, time_point, sample,
                                group_number, `repeat`) %>% 
  summarize(read_count = n())

# Analyze taxonomy of each CRISPR array (group) at the genus level
genus <- total %>% group_by(individual, time_point, sample, 
                                  group_number, `repeat`, phylum, genus) %>% 
  summarize(genus_count = n())
genus <- left_join(genus, group_sum)
genus <- genus %>% mutate(genus_content = genus_count/read_count)
genus <- genus %>% mutate(sample_group = paste(sample, group_number, sep = '_'))

genus_represent=data.frame()
unique(genus$sample_group)

for (i in unique(genus$sample_group)) {
  print(i)
  temp_df <- genus[grep(i, genus$sample_group), ]
  temp_df <- temp_df %>% arrange(desc(genus_content))
  genus_represent <- rbind(genus_represent, temp_df[1, ])
  temp_df = data.frame()
}

# Analyze taxonomy of each CRISPR array (group) at the phylum level
phylum <- total %>% group_by(individual, time_point, sample, 
                             group_number, `repeat`, phylum) %>% 
  summarize(phylum_count = n())
phylum <- left_join(phylum, group_sum)
phylum <- phylum %>% mutate(phylum_content = phylum_count/read_count)
phylum <- phylum %>% mutate(sample_group = paste(sample, group_number, sep = '_'))

phylum_represent=data.frame()
for (i in unique(phylum$sample_group)){
  print(i)
  temp_df <- phylum[phylum$sample_group == i, ]
  phylum_represent <- rbind(phylum_represent, temp_df[1, ])
  temp_df <- data.frame()
}

################################################################################
genus_plot <- genus_represent %>% group_by(individual, time_point, genus) %>% 
  summarise(genus_sum = n())
genus_plot %>% filter(genus != 'unclassified') %>% ggplot(aes(time_point, genus_sum)) +
  facet_grid(.~individual) + 
  geom_col(aes(fill = genus)) + 
  theme(legend.position = 'None')

stat(total$classified)

################################################################################
total_0.15 <- read.csv("../intermediates/17.kraken_taxonomy/confidence_0.15/kraken_total_output.csv", 
                      sep = ',', header = F)
total_0.15 <- total_0.15 %>% separate(V2, c("individual", "time_point", "sample", "group", 
                                  "group_number","repeat", "SRR", "read"))
colnames(total_0.15)[c(1, 10, 11, 12, 13)] <- c('classified', 'taxonomy', 'taxid', 
                                           'read_length', 'kmer_mapping')
evenness_0.15 <- total_0.15 %>% filter(classified == 'C') %>% 
  group_by(individual, time_point, sample, group_number, taxonomy) %>% 
  summarise(count = n())
evenness_0.15 <- evenness_0.15 %>% group_by(individual, time_point, sample, group_number) %>% 
  summarise(count = n())

total_0.6 <- read.csv("../intermediates/17.kraken_taxonomy/confidence_0.6/kraken_total_output.csv", 
                       sep = ',', header = F)
total_0.6 <- total_0.6 %>% separate(V2, c("individual", "time_point", "sample", "group", 
                                            "group_number","repeat", "SRR", "read"))
colnames(total_0.6)[c(1, 10, 11, 12, 13)] <- c('classified', 'taxonomy', 'taxid', 
                                                'read_length', 'kmer_mapping')
evenness_0.6 <- total_0.6 %>% filter(classified == 'C') %>% 
  group_by(individual, time_point, sample, group_number, taxonomy) %>% 
  summarise(count = n())
evenness_0.6 <- evenness_0.6 %>% group_by(individual, time_point, sample, group_number) %>% 
  summarise(count = n())
################################################################################

unclassified_detection <- as.data.frame(matrix(ncol = 7))
colnames(unclassified_detection) <- c('kingdom', 'phylum', 'class', 'order', 
                                      'family', 'genus', 'species')
for (i in 1:nrow(total[1:100, 13:19])){
  unclassified_detection[i,] <- str_detect(total[i, 13:19], "unclassified")
  print(i)
}
################################################################################


# Run code from here
library(stringr)
library(tidyr)
library(dplyr)
getwd()
#-------------------------------------------------------------------------------
# Determine the classification uniqueness of different confidence score (0.15 and 0.6) with Kraken.

# confidence score 0.15
total_0.15 <- read.csv('../intermediates/17.kraken_taxonomy/confidence_0.15/kraken_total_output.csv', 
                       sep = ',', header = F)
total_0.15 <- total_0.15 %>% separate(V2, c("individual", "time_point", "sample", "group", 
                                            "group_number","repeat", "SRR", "read"))
colnames(total_0.15)[c(1, 10, 11, 12, 13)] <- c('classified', 'taxonomy', 'taxid', 
                                               'read_length', 'kmer_mapping')
total_0.15 <- total_0.15[, -5]
taxid_0.15 <- total_0.15$taxid
write.csv(taxid_0.15, file = '../intermediates/17.kraken_taxonomy/confidence_0.15/taxid.txt', 
          quote = F,row.names = F, col.names = F, sep = '\t')
taxonomy_0.15 <- read.table('../intermediates/17.kraken_taxonomy/confidence_0.15/taxonomy.txt', 
                            header = F, sep = '\t', fill = NA)
taxonomy_0.15 <- taxonomy_0.15[-1,]
taxonomy_0.15 <- taxonomy_0.15[,-2]
taxonomy_0.15[taxonomy_0.15 == ''] <- 'unclassified'
colnames(taxonomy_0.15) <- c("taxid", "kingdom", "phylum", "class", "order", 
                             "family", "genus", "species")
for (i in 2:ncol(taxonomy_0.15)){
  taxonomy_0.15[, i][str_detect(taxonomy_0.15[, i], 'unclassified')] <- "unclassified"
  print(i)
}
total_0.15 <- cbind(total_0.15, taxonomy_0.15[, 2:8])
write.table(total_0.15, file = '../intermediates/17.kraken_taxonomy/confidence_0.15/kraken_total_summary.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
genus_0.15 <- total_0.15 %>% group_by(individual, time_point, sample, group_number,
                                      genus) %>% 
  summarise(genus_count = n())
genus_uniqueness_0.15 <- genus_0.15 %>% group_by(individual, time_point, sample, 
                                                 group_number) %>% 
  summarise(genus_uniqueness = n())


# confidence score 0.6
total_0.6 <- read.csv('../intermediates/17.kraken_taxonomy/confidence_0.6/kraken_total_output.csv', 
                       sep = ',', header = F)
total_0.6 <- total_0.6 %>% separate(V2, c("individual", "time_point", "sample", "group", 
                                            "group_number","repeat", "SRR", "read"))
colnames(total_0.6)[c(1, 10, 11, 12, 13)] <- c('classified', 'taxonomy', 'taxid', 
                                                'read_length', 'kmer_mapping')
total_0.6 <- total_0.6[, -5]
taxid_0.6 <- total_0.6$taxid
write.csv(taxid_0.6, file = '../intermediates/17.kraken_taxonomy/confidence_0.6/taxid.txt', 
          quote = F,row.names = F, col.names = F, sep = '\t')
taxonomy_0.6 <- read.csv('../intermediates/17.kraken_taxonomy/confidence_0.6/taxonomy.txt', 
                            header = F, sep = '\t', fill = NA) # due to too many blanks in the taxonomy.txt, \t was manually added at row1 to allow the read.csv to read the table into 9 columns.
taxonomy_0.6 <- taxonomy_0.6[-1,]
taxonomy_0.6 <- taxonomy_0.6[,-2]
taxonomy_0.6[taxonomy_0.6 == ''] <- 'unclassified'
colnames(taxonomy_0.6) <- c("taxid", "kingdom", "phylum", "class", "order", 
                             "family", "genus", "species")
for (i in 2:ncol(taxonomy_0.6)){
  taxonomy_0.6[, i][str_detect(taxonomy_0.6[, i], 'unclassified')] <- "unclassified"
  print(i)
}
total_0.6 <- cbind(total_0.6, taxonomy_0.6[, 2:8])
write.table(total_0.6, file = '../intermediates/17.kraken_taxonomy/confidence_0.6/kraken_total_summary.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
genus_0.6 <- total_0.6 %>% group_by(individual, time_point, sample, group_number,
                                      genus) %>% 
  summarise(genus_count = n())
genus_uniqueness_0.6 <- genus_0.6 %>% group_by(individual, time_point, sample, 
                                                 group_number) %>% 
  summarise(genus_uniqueness = n())