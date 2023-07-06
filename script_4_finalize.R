library(tidyverse)
library(microseq)

args <- commandArgs(trailingOnly = T)
contigs_folder <- args[1]

# contigs_folder <- "iter1/contigs"

contigs_files <- list.files(contigs_folder, pattern = "_contigs.fasta", full.names = T)
estimated_16S.tbl <- NULL
for(i in 1:length(contigs_files)){
  seed_name <- str_remove(basename(contigs_files[i]), "_contigs.fasta")
  estimated_16S.tbl <- readFasta(contigs_files[i]) %>% 
    mutate(length = str_length(Sequence)) %>% 
    arrange(desc(length)) %>% 
    slice(1) %>% 
    select(-length) %>% 
    mutate(Header = str_c(seed_name, " ", Header)) %>% 
    bind_rows(estimated_16S.tbl)
}
writeFasta(estimated_16S.tbl, out.file = file.path(contigs_folder, "estimated_16S.fasta"))
