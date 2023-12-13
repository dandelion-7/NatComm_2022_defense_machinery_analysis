~/crisprome/human_temporal_metagenome_wangjun_2022natComm/scripts

1.fastqc_1.sh
This script is for primarily perform quality control of raw metagenomic sequencing data.

2.fastp.sh
Following script 1, this is for filtering low-quality reads, and removing adaptors.

3.fastqc_2.sh
Following script 2, this is for re-check the quality of processed sequencing data.

4.bowtie2_dehost.sh
Following script 2, this is for removing the host(human)-derived reads in the processed data.

5.crass.sh
Following script 4, crass will detect all the possible CRISPR arrays in the metagenomic reads and extract the correspoding reads.

6.metaspades.sh
Following script 4, metaspades assembled each individual's sequencing files at each time point into contigs.

7.crisprtools_extract.sh
Following script 5, use crisprtools to extract all identified repeats, spacers, and flankers.

8.cdhit.sh
Following script 7, this is for clustering each individual's identified CRISPR arrays to define the stable CRISPR arrays.

9.bowtie2.sh
Following script 10, use bowtie2 to build reference with each one's contigs and align sequencing data to the reference to determine the abundancies of all the assembled contigs, which is the prerequisite for binning.

10.megahit.sh
Following script 4, this script merges each individual's data at ten time points and assemble into contigs.
metaspades can't assemble an individual's datasets from ten time points due to the limitation of RAW. So megahit is used instead of 6.metaspades.sh.

11.clustered repeats analysis
Following script 7 and 8, the clstuering results of recovered CRISPRs are summarized (11-1). Kmers of the repeats are calculated (11-2), and with the kmers of repeats, each individuals' repeats kmer profile are used to calculate the Bray-Curtis distances and NMDS to show the clustering/distances of each one's repeats.

12.metabat2.sh
Following script 9 and 10, with the contigs assembled from megahits, and sorted bam files from bowtie2 through mapping reads to the contigs, the contigs are binned with metabat2.

13.repeats blastn
Following script 7, with the repeats database constructed from 5 published datasets, the recovered repeats are blasted against the database to verify if they are real CRISPR arrays (already discovered in other studies).

14.checkm
Following script 12, CheckM is used to evalutae the integrity and abundance of each bin.

15.mmseqs_taxonomy
Following script 10, MMseqs2 taxonomy is used to assign taxonomy to the assembled contigs.

16.bowtie2_repeats_contig
Following script 7/10, repeat sequences are aligned to assembled contigs with bowtie2. These contigs have been assiged taxonomic information with MMseqs2, so that the taxonomic source of these CRISPR arrays can be identified.

17.kraken_taxonomy
In parallele with script 15, and following script 5, kraken2 is directly used to assign taxonomies to the reads corresponding to CRISPRs recovered by crass in script 5.


