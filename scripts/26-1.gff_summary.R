getwd()
setwd('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts/')

# summarize all the CDSs predicted with Prodigal.
total_CDS_info <- read.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/25.defensefinder/prodigal/total_CDS.txt', 
                   sep = ';', header = F)
colnames(total_CDS_info) <- c('CDS', 'start', 'end', 'strand', 'ID', 'partial', 
                              'start_type', 'rbs_motif', 'rbs_spacer', 'gc_cont')
total_CDS_info <- total_CDS_info %>% mutate(attributes = paste(ID, partial, start_type, rbs_motif, rbs_spacer, gc_cont, 
                                             sep = ';'))
total_CDS_info <- total_CDS_info[, -c(5:10)]

CDSs <- as.data.frame(total_CDS_info$CDS)
CDSs <- CDSs %>% separate(`total_CDS_info$CDS`, sep = '_', 
                          into = c('individual', 'k141', 'contig', 'flag', 'multi', 'len', 'CDS_code'))
CDSs <- CDSs %>% mutate(contig_id = paste(k141, contig, sep = '_'))
total_CDS_info <- cbind(total_CDS_info, CDSs[, c(1, 8, 6, 7)])

# summarize all the anti-phage systems predicted with DefenseFinder
total_antiphage_systems <- read.delim('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/25.defensefinder/defensefinder/merged_defensefinder_systems.tsv', 
                                      sep = '\t', header = T)
system_size <- total_antiphage_systems[, c('sys_id', 'genes_count')]

# summarize all the anti-phage genes predicted with DefenseFinder
total_antiphage_genes <- read.delim('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/25.defensefinder/defensefinder/merged_defensefinder_genes.tsv', 
                                    sep = '\t', header = T)
colnames(total_antiphage_genes)[2] <- 'CDS'
total_antiphage_genes <- left_join(total_antiphage_genes, system_size, by = 'sys_id')

# merge lists of anti-phage genes and total CDSs to generate annotation.
annotation_list <- left_join(total_CDS_info, total_antiphage_genes, by = 'CDS')
annotation_list <- annotation_list[, -c(10, 29, 30)]
annotation_list$gene_name[is.na(annotation_list$gene_name)] <- 'non-antiphage-gene'
annotation_list$gene_name <- paste("gene_name=", annotation_list$gene_name, sep = '')
annotation_list$CDS <- str_replace_all(annotation_list$CDS, pattern = '=', replacement = '^')
# all the "=" in the CDS id are replaced with ^ because in the attributes of gff, no extra "=" are allowed.
annotation_list$gene_name <- paste(annotation_list$gene_name, annotation_list$CDS, sep = '+')
annotation_list$genes_count[is.na(annotation_list$genes_count)] <- 0
annotation_list$system_size <- paste("system_size=", annotation_list$genes_count, sep = '')
annotation_list$attributes <- paste(annotation_list$attributes, annotation_list$gene_name, annotation_list$system_size, sep = ";")
annotation_list[is.na(annotation_list)] <- '.'
annotation_list$source <- 'Prodigal+DefenseFinder'
annotation_list$phase <- 0
annotation_list$type <- 'CDS'
annotation_list %>% write.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/26-1.total_CDSs_annotation_list.txt', 
                                sep = '\t', quote = F, row.names = F, col.names = T)

# generate GFF file.
gff <- annotation_list[, c(6, 7, 8, 9, 1, 28, 30, 2, 3, 22, 4, 29, 5)]
total_genes_gff <- annotation_list[, c('contig_id', 'source', 'type', 'start', 
                                       'end', 'hit_i_eval', 'strand', 'phase', 
                                       'attributes', 'individual')]
total_genes_gff %>% write.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/gff/26-1.total_CDSs_gff.txt', 
                                sep = '\t', quote = F, row.names = F, col.names = F)
gff_1 <- filter(total_genes_gff, individual == 1)[, -10]
gff_1 %>% write.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/gff/26-1.1.gff', 
                      sep = '\t', quote = F, row.names = F, col.names = F)

gff_2 <- filter(total_genes_gff, individual == 2)[, -10]
gff_2 %>% write.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/gff/26-1.2.gff', 
                      sep = '\t', quote = F, row.names = F, col.names = F)

gff_3 <- filter(total_genes_gff, individual == 3)[, -10]
gff_3 %>% write.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/gff/26-1.3.gff', 
                      sep = '\t', quote = F, row.names = F, col.names = F)

gff_4 <- filter(total_genes_gff, individual == 4)[, -10]
gff_4 %>% write.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/gff/26-1.4.gff', 
                      sep = '\t', quote = F, row.names = F, col.names = F)

gff_5 <- filter(total_genes_gff, individual == 5)[, -10]
gff_5 %>% write.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/gff/26-1.5.gff', 
                      sep = '\t', quote = F, row.names = F, col.names = F)

gff_6 <- filter(total_genes_gff, individual == 6)[, -10]
gff_6 %>% write.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/gff/26-1.6.gff', 
                      sep = '\t', quote = F, row.names = F, col.names = F)

gff_7 <- filter(total_genes_gff, individual == 7)[, -10]
gff_7 %>% write.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/gff/26-1.7.gff', 
                      sep = '\t', quote = F, row.names = F, col.names = F)

gff_8 <- filter(total_genes_gff, individual == 8)[, -10]
gff_8 %>% write.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/gff/26-1.8.gff', 
                      sep = '\t', quote = F, row.names = F, col.names = F)

gff_9 <- filter(total_genes_gff, individual == 9)[, -10]
gff_9 %>% write.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/gff/26-1.9.gff', 
                      sep = '\t', quote = F, row.names = F, col.names = F)

gff_10 <- filter(total_genes_gff, individual == 10)[, -10]
gff_10 %>% write.table('~/crisprome/human_temporal_metagenome_wangjun_2022natComm/intermediates/26.antivirus_genes_analysis/gff/26-1.10.gff', 
                      sep = '\t', quote = F, row.names = F, col.names = F)
