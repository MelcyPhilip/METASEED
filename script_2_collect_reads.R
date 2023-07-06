library(tidyverse)
library(microseq)


args <- commandArgs(trailingOnly = T)
prefix <-args[1]
seed_abundance_file <- args[2]
r1_fastq_file <- args[3]
r2_fastq_file <- args[4]
aln_folder <-args[5]
seed_reads_folder <-args[6]


#############
### Settings
###
seqkit_exe <- "/mnt/users/larssn/bin/seqkit"
max_n <- 3


##########################
### The seed abundances
###
seed.tbl <- suppressMessages(read_delim(seed_abundance_file, delim = "\t")) %>%
  select(1, all_of(prefix))
colnames(seed.tbl) <- c("seed_id", "abundance")


###########################
### The blast alignments
###
r1_blast_file <- file.path(aln_folder, str_c(prefix, "_R1.txt"))
r2_blast_file <- file.path(aln_folder, str_c(prefix, "_R2.txt"))
cat("Sample", prefix, "\n")
b1.tbl <- suppressMessages(read_delim(r1_blast_file, delim = "\t",
                                      col_names = c("read_id", "seed_id", "pident", "gapopen",
                                                    "qlen", "qstart", "qend",
                                                    "slen", "sstart", "send", "bitscore"))) %>% 
  mutate(type = "R1")
b2.tbl <- suppressMessages(read_delim(r2_blast_file, delim = "\t",
                                      col_names = c("read_id", "seed_id", "pident", "gapopen",
                                                    "qlen", "qstart", "qend",
                                                    "slen", "sstart", "send", "bitscore"))) %>% 
  mutate(type = "R2")
b.tbl <- bind_rows(b1.tbl, b2.tbl) %>% 
  mutate(readpair_id = str_remove(read_id, "/[12]$")) %>%  # read 1 and 2 may have /1 and /2 at the end of the name
  filter(gapopen == 0) %>% 
  mutate(strand = if_else(sstart < send, '+', '-')) %>% 
  mutate(match = if_else(qstart == 1 & qend == qlen, 'inside', 'local')) %>%
  mutate(match = if_else((qstart == 1 & qend < qlen) & ((send == slen & strand == '+') | (send == 1 & strand == '-')), 'start_end', match)) %>%
  mutate(match = if_else((qstart > 1 & qend == qlen) & ((sstart == 1 & strand == '+') | (sstart == slen & strand == '-')), 'end_start', match)) %>% 
  filter(match != "local") %>% 
  distinct(read_id, seed_id, .keep_all = T) %>%        # only one combination of read and seed
  left_join(seed.tbl, by = "seed_id") %>%              # include abundance column
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

rm(b1.tbl, b2.tbl)
u.seeds <- unique(b.tbl$seed_id)
cat("   found", nrow(b.tbl), "alignments matching", length(u.seeds), "unique seeds...\n")

cat("   collecting R1 seed-reads...\n")
b.tbl %>% 
  mutate(read_id = str_replace(read_id, "/2$", "/1")) %>% 
  distinct(read_id) %>% 
  write_delim(delim = "\t", file = file.path(aln_folder, str_c(prefix, "_tags.txt")), col_names = F)
cmd <- paste(seqkit_exe,
             "grep -f", file.path(aln_folder, str_c(prefix, "_tags.txt")),
             r1_fastq_file, ">",
             file.path(aln_folder, str_c(prefix, "_R1.fq")))
system(cmd)
fq1 <- readFastq(file.path(aln_folder, str_c(prefix, "_R1.fq"))) %>% 
  mutate(Header = word(Header, 1))
cat("   collected", nrow(fq1), "reads from the R1 file...\n")
for(i in 1:length(u.seeds)){
  s.tbl <- b.tbl %>% 
    filter(seed_id == u.seeds[i]) %>% 
    mutate(read_id = str_replace(read_id, "/2$", "/1"))
  fq1 %>% 
    filter(Header %in% s.tbl$read_id)%>% 
    writeFastq(file.path(seed_reads_folder, str_c(prefix, "_", u.seeds[i], "_R1.fq")))
}
ok <- file.remove(file.path(aln_folder, str_c(prefix, "_R1.fq")))

cat("   collecting R2 seed-reads...\n")
b.tbl %>% 
  mutate(read_id = str_replace(read_id, "/1$", "/2")) %>% 
  distinct(read_id) %>% 
  write_delim(delim = "\t", file = file.path(aln_folder, str_c(prefix, "_tags.txt")), col_names = F)
cmd <- paste(seqkit_exe,
             "grep -f", file.path(aln_folder, str_c(prefix, "_tags.txt")),
             r2_fastq_file, ">",
             file.path(aln_folder, str_c(prefix, "_R2.fq")))
system(cmd)
fq2 <- readFastq(file.path(aln_folder, str_c(prefix, "_R2.fq"))) %>% 
  mutate(Header = word(Header, 1))
cat("   collected", nrow(fq2), "reads from the R2 file...\n")
for(i in 1:length(u.seeds)){
  s.tbl <- b.tbl %>% 
    filter(seed_id == u.seeds[i]) %>% 
    mutate(read_id = str_replace(read_id, "/1$", "/2"))
  fq2 %>% 
    filter(Header %in% s.tbl$read_id) %>% 
    writeFastq(file.path(seed_reads_folder, str_c(prefix, "_", u.seeds[i], "_R2.fq")))
}
ok <- file.remove(file.path(aln_folder, str_c(prefix, "_R2.fq")))
ok <- file.remove(file.path(aln_folder, str_c(prefix, "_tags.txt")))
cat("done!\n")
