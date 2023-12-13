# This script is for summarizing the blast results of recovered repeats to repeats from data. Not finished yet.
#----------------------------------------------------------------------------------------------------------------------------------
getwd()
blastn_results <- read.table('../intermediates/13.repeats_blast_on_database/blastn_results/total_blastn_database.txt', 
                             header = TRUE, sep = '\t')
blastn_results <- blastn_results %>% mutate(subject_coverage = nident / slen)
blastn_results <- blastn_results %>% mutate(query_subject_length = qlen - slen)
