#!/bin/bash
sample_id=$1
R1_file=$2
R2_file=$3
dbase=$4
tmp_folder=$5
seed_abundance_file=$6
max_n=$7
threads=$8


echo "*** Mapping $R1_file reads to seeds with bowtie2, and collecting matching reads"
$bowtie2_app bowtie2 --threads $threads --local -x $dbase -U $R1_file \
 | $samtools_app samtools view --threads $threads -b - \
 | $samtools_app samtools fasta -0 $tmp_folder/$sample_id\_R1.fasta -F 4 -
echo "*** Mapping $R2_file reads to seeds with bowtie2, and collecting matching reads"
$bowtie2_app bowtie2 --threads $threads --local -x $dbase -U $R2_file \
 | $samtools_app samtools view --threads $threads -b - \
 | $samtools_app samtools fasta -0 $tmp_folder/$sample_id\_R2.fasta -F 4 -


echo "*** Aligning collected R1-reads to seeds with BLAST..."
$blast_app blastn -db $dbase -query $tmp_folder/$sample_id\_R1.fasta -perc_identity 99 -num_threads $threads \
  -outfmt '6 qseqid sseqid pident gapopen qlen qstart qend slen sstart send bitscore' \
  -penalty -1 -reward 1 -gapopen 2 -gapextend 1 -out $tmp_folder/$sample_id\_R1.txt
echo "*** Aligning collected R1-reads to seeds with BLAST..."
$blast_app blastn -db $dbase -query $tmp_folder/$sample_id\_R2.fasta -perc_identity 99 -num_threads $threads \
  -outfmt '6 qseqid sseqid pident gapopen qlen qstart qend slen sstart send bitscore' \
  -penalty -1 -reward 1 -gapopen 2 -gapextend 1 -out $tmp_folder/$sample_id\_R2.txt


echo "*** Reading BLAST results and finding the matching read-pairs"
Rscript -e "suppressMessages(library(tidyverse))
args <- commandArgs(trailingOnly = T)
seed_abundance_file <- args[1]
tmp_folder <- args[2]
sample_id <- args[3]
max_n <- args[4]
seed.tbl <- suppressMessages(read_delim(seed_abundance_file, delim = '\t')) %>%
  select(1, all_of(sample_id))
colnames(seed.tbl) <- c('seed_id', 'abundance')
r1_blast_file <- file.path(tmp_folder, str_c(sample_id, '_R1.txt'))
r2_blast_file <- file.path(tmp_folder, str_c(sample_id, '_R2.txt'))
b1.tbl <- suppressMessages(read_delim(r1_blast_file, delim = '\t',
                                      col_names = c('read_id', 'seed_id', 'pident', 'gapopen',
                                                    'qlen', 'qstart', 'qend',
                                                    'slen', 'sstart', 'send', 'bitscore'))) %>%
  mutate(type = 'R1')
b2.tbl <- suppressMessages(read_delim(r2_blast_file, delim = '\t',
                                      col_names = c('read_id', 'seed_id', 'pident', 'gapopen',
                                                    'qlen', 'qstart', 'qend',
                                                    'slen', 'sstart', 'send', 'bitscore'))) %>%
  mutate(type = 'R2')
b.tbl <- bind_rows(b1.tbl, b2.tbl) %>%
  mutate(readpair_id = str_remove(read_id, '/[12]$')) %>%  # read 1 and 2 may have /1 and /2 at the end of the name
  filter(gapopen == 0) %>%
  mutate(strand = if_else(sstart < send, '+', '-')) %>%
  mutate(match = if_else(qstart == 1 & qend == qlen, 'inside', 'local')) %>%
  mutate(match = if_else((qstart == 1 & qend < qlen) & ((send == slen & strand == '+') | (send == 1 & strand == '-')), 'start_end', match)) %>%
  mutate(match = if_else((qstart > 1 & qend == qlen) & ((sstart == 1 & strand == '+') | (sstart == slen & strand == '-')), 'end_start', match)) %>%
  filter(match != 'local') %>%
  distinct(read_id, seed_id, .keep_all = T) %>%        # only one combination of read and seed
  left_join(seed.tbl, by = 'seed_id') %>%              # include abundance column
  arrange(desc(abundance)) %>%                         # arrange by abundance, descending order
  group_by(read_id) %>%                                # group by read_id
  mutate(max_score = max(bitscore)) %>%                # compute max_score for each group
  filter(bitscore == max_score) %>%                    # keep only alignments with max_score
  mutate(n_hits = length(unique(seed_id))) %>%         # then count how many different seeds a read matches
  filter(n_hits <= max_n) %>%                          # keep only those with max_n hits or less
  slice(1) %>%                                         # keep only the largest abundance within group
  ungroup() %>%
  group_by(readpair_id) %>%                            # in case the two reads in a pair disagree on seed
  slice(1) %>%                                         # again we keep the seed with largest abundance
  ungroup() %>%
  arrange(desc(n_hits), readpair_id, desc(abundance)) %>%
  select(readpair_id, read_id, seed_id, n_hits, abundance)
cat('Found', nrow(b.tbl), 'alignments matching', length(unique(b.tbl\$seed_id)), 'unique seeds...\n')
b.tbl %>%
  mutate(read_id = str_replace(read_id, '/2$', '/1')) %>%
  distinct(read_id) %>%
  write_delim(delim = '\t', file = file.path(tmp_folder, str_c(sample_id, '_R1_tags.txt')), col_names = F)
b.tbl %>%
  mutate(read_id = str_replace(read_id, '/1$', '/2')) %>%
  distinct(read_id) %>%
  write_delim(delim = '\t', file = file.path(tmp_folder, str_c(sample_id, '_R2_tags.txt')), col_names = F)
write_delim(b.tbl, delim = '\t', file = file.path(tmp_folder, str_c(sample_id, '_b.tbl')))" $seed_abundance_file $tmp_folder $sample_id $max_n


echo "*** Collecting the actual read-pairs"
$bbmap_app filterbyname.sh in=$R1_file out=$tmp_folder/$sample_id\_R1.fq names=$tmp_folder/$sample_id\_R1_tags.txt include=t
$bbmap_app filterbyname.sh in=$R2_file out=$tmp_folder/$sample_id\_R2.fq names=$tmp_folder/$sample_id\_R2_tags.txt include=t


echo "*** Distributing reads by seeds"
Rscript -e "suppressMessages(library(tidyverse))
suppressMessages(library(microseq))
args <- commandArgs(trailingOnly = T)
tmp_folder <- args[1]
sample_id <- args[2]
b.tbl <- read_delim(file.path(tmp_folder, str_c(sample_id, '_b.tbl')))
fq1 <- readFastq(file.path(tmp_folder, str_c(sample_id, '_R1.fq'))) %>%
  mutate(Header = word(Header, 1))
fq2 <- readFastq(file.path(tmp_folder, str_c(sample_id, '_R2.fq'))) %>%
  mutate(Header = word(Header, 1))
u.seeds <- unique(b.tbl\$seed_id)
for(i in 1:length(u.seeds)){
  s.tbl <- b.tbl %>%
    filter(seed_id == u.seeds[i]) %>%
    mutate(read_id = str_replace(read_id, '/2$', '/1'))
  fq1 %>%
    filter(Header %in% s.tbl\$read_id)%>%
    writeFastq(file.path(tmp_folder, str_c(sample_id, '_', u.seeds[i], '_R1.fq')))
  s.tbl <- b.tbl %>%
    filter(seed_id == u.seeds[i]) %>%
    mutate(read_id = str_replace(read_id, '/1$', '/2'))
  fq2 %>%
    filter(Header %in% s.tbl\$read_id) %>%
    writeFastq(file.path(tmp_folder, str_c(sample_id, '_', u.seeds[i], '_R2.fq')))
}" $tmp_folder $sample_id



