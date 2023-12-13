getwd()
setwd('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/')

#-------------------------------------------------------------------------------
# analyse results of blasting spacers onto predicted virus contigs.
virus <- read.table('../intermediates/21.spacers_blast_on_targets/spacers_blast_on_mixed_virus_plasmid/all_smaples_virus_contigs.txt', 
                    sep = '\t', header = F)
plasmid <- read.table('../intermediates/21.spacers_blast_on_targets/spacers_blast_on_mixed_virus_plasmid/all_smaples_plasmid_contigs.txt', 
                      sep = '\t', header = F)
colnames(virus) <- c('qseqid', 'sseqid', 'pident', 'nident', 'qlen', 'slen',
                     'length', 'mismatch', 'positive', 'ppos', 'gapopen',
                     'gaps', 'qstart', 'qend', 'sstart', 'send', 'evalue', 
                     'bitscore', 'qcovs', 'qcovhsp', 'qcovus', 'qseq', 'sstrand',
                     'frames')
colnames(plasmid) <- c('qseqid', 'sseqid', 'pident', 'nident', 'qlen', 'slen',
                     'length', 'mismatch', 'positive', 'ppos', 'gapopen',
                     'gaps', 'qstart', 'qend', 'sstart', 'send', 'evalue', 
                     'bitscore', 'qcovs', 'qcovhsp', 'qcovus', 'qseq', 'sstrand',
                     'frames')
virus$protospacer_origin <- 'virus'
plasmid$protospacer_origin <- 'plasmid'
total <- rbind(virus, plasmid)

query <- as.data.frame(total$qseqid)
query <- query %>% separate(`total$qseqid`, sep = '_',
                            into = c('individual', 'time_point', 'sample', 
                                     'CRISPR_group', 'Cov', 'coverage'))
query <- query %>% separate(CRISPR_group, sep = 'SP', into = c('CRISPR_group', 'spacer_code'))
query <- query[, -6]

ref <- as.data.frame(total$sseqid)
ref <- ref %>% separate(`total$sseqid`, sep = '_', into = c('contig_individual'))

total <- cbind(query, ref, total)
total <- total %>% filter(qcovs == 100 & pident > 95)

spacer_protospacer <- total %>% group_by(individual, contig_individual) %>% 
  summarise(pair_count = n())
spacer_sum <- spacer_protospacer %>% group_by(individual) %>% 
  summarise(total_pair = sum(pair_count))
spacer_protospacer <- left_join(spacer_protospacer, spacer_sum)
spacer_protospacer <- spacer_protospacer %>% mutate(
  pair_proportion = pair_count/total_pair
)

spacer_protospacer$individual <- as.numeric(spacer_protospacer$individual)
spacer_protospacer$contig_individual <- factor(spacer_protospacer$contig_individual, 
                                               levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
color_10 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
              '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')
spacer_protospacer %>% ggplot(aes(x = individual, y = pair_proportion)) + 
  geom_col(aes(fill = contig_individual), 
           position = position_dodge()) + 
  scale_fill_manual(values = color_10) + 
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

possible_phage_CRISPR <- total %>% group_by(individual, time_point, sample, CRISPR_group, sseqid) %>% 
  summarise(repetitive_targeting = n())
possible_phage_CRISPR %>% ggplot(aes(x = repetitive_targeting)) + 
  geom_histogram() + 
  scale_x_continuous(limits = c(0, 10))

total <- total %>% mutate(targeting_pair = paste(individual, time_point, CRISPR_group, sseqid, sep = '-'))
total <- left_join(total, possible_phage_CRISPR[, 6:7])
total %>% ggplot(aes(x = repetitive_targeting)) + 
  geom_histogram() + 
  facet_grid(protospacer_origin ~ .)
total_filtered <- total %>% filter(protospacer_origin == 'virus')

spacer_protospacer <- total_filtered %>% group_by(individual, contig_individual) %>% 
  summarise(pair_count = n())
spacer_sum <- spacer_protospacer %>% group_by(individual) %>% 
  summarise(total_pair = sum(pair_count))

# summarize the total number of virus contigs from each individual
contig_sum <- read.table('../intermediates/21.spacers_blast_on_targets/spacers_blast_on_mixed_virus_plasmid/blastdb/all_individuals/individual_virus_contigs_count.txt', 
                         sep = ' ', header = F)[,4:5]
colnames(contig_sum) <- c('contig_count', 'contig_individual')
contig_sum$contig_individual <- as.character(contig_sum$contig_individual)

# summarize the total number of spacers from each individual
spacer_sum <- read.table('../intermediates/7.crisprtools_extract/all_spacers_list.txt', sep = '_', header = F)
spacer_sum$V1 <- str_replace(spacer_sum$V1, '>', '')
spacer_sum <- spacer_sum[, -5]
colnames(spacer_sum) <- c('individual', 'time_point', 'sample', 'CRISPR_group', 'coverage')
spacer_sum <- spacer_sum %>% group_by(individual) %>% 
  summarize(spacer_count = n())
colnames(spacer_sum) <- c('individual', 'spacer_count')

spacer_protospacer <- left_join(spacer_protospacer, spacer_sum, by = 'individual')
spacer_protospacer <- left_join(spacer_protospacer, contig_sum, by = 'contig_individual')
spacer_protospacer <- spacer_protospacer %>% mutate(
  pair_proportion = pair_count/spacer_count/contig_count
)
spacer_protospacer$individual <- factor(spacer_protospacer$individual, 
                                        levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
spacer_protospacer$contig_individual <- factor(spacer_protospacer$contig_individual, 
                                               levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
color_10 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99',
              '#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')
spacer_protospacer %>% ggplot(aes(x = contig_individual, y = pair_proportion)) + 
  geom_col(aes(fill = individual), 
           position = position_dodge()) + 
  scale_fill_manual(values = color_10) +  
  scale_x_discrete(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
spacer_protospacer %>% ggplot(aes(x = individual, y = pair_proportion)) + 
  geom_col(aes(fill = contig_individual), 
           position = position_dodge(), 
           width = 0.5) + 
  scale_fill_manual(values = color_10) +  
  scale_x_discrete(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)) + 
  ylab('#targeting pair / (#individual virus contig ✖ #individual spacer)') + 
  theme(legend.position = 'top') + 
  guides(fill = guide_legend(nrow = 1)) + 
  xlab('individual belonging of spacers') + 
  labs(fill = 'individual belonging of virus contigs')
# ~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/21.spacers_blast_on_targets/results_analysis/figures/21-4-1.individual_spacer_targeting_preference

#-------------------------------------------------------------------------------
# analyse results of blast spacers onto public virome databases.
virus_seq <- read.table('../intermediates/21.spacers_blast_on_targets/spacers_blast_on_database/all_samples_virus_sequences.txt', 
                        sep = '\t', header = F)
huvirdb <- read.table('../intermediates/21.spacers_blast_on_targets/spacers_blast_on_database/all_samples_HuVirDB.txt', 
                      sep = '\t', header = F)
colnames(virus_seq) <- c('qseqid', 'sseqid', 'pident', 'nident', 'qlen', 'slen',
                         'length', 'mismatch', 'positive', 'ppos', 'gapopen',
                         'gaps', 'qstart', 'qend', 'sstart', 'send', 'evalue', 
                         'bitscore', 'qcovs', 'qcovhsp', 'qcovus', 'qseq', 'sstrand',
                         'frames')
colnames(huvirdb) <- c('qseqid', 'sseqid', 'pident', 'nident', 'qlen', 'slen',
                       'length', 'mismatch', 'positive', 'ppos', 'gapopen',
                       'gaps', 'qstart', 'qend', 'sstart', 'send', 'evalue', 
                       'bitscore', 'qcovs', 'qcovhsp', 'qcovus', 'qseq', 'sstrand',
                       'frames')

virus_seq$ref <- 'mixed_virome_db'
huvirdb$ref <- 'huvirdb'
total <- rbind(virus_seq, huvirdb)
total <- total %>% filter(qcovs == 100 & pident > 95)
db_diff_1 <- total %>% group_by(qseqid, ref) %>%
 summarise(targeting = 1, target_count = n())
db_diff_2 <- db_diff_1 %>% group_by(qseqid) %>%
  summarize(targeted_db = n_distinct(ref))
summary(db_diff_2$targeted_db)
