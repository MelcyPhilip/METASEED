#!/bin/bash
seed_sequence_file=$1
tmp_folder=$2
out_folder=$3
threads=$4
dbase=$5


echo "*** Merging reads for the same seed from different samples"
echo "*** into a single pair of fastq files"
Rscript -e "suppressMessages(library(microseq))
args <- commandArgs(trailingOnly = T)
seed_sequence_file <- args[1]
tmp_folder <- args[2]
seeds.tbl <- readFasta(seed_sequence_file) %>%
  mutate(Header = word(Header, 1)) %>%
  mutate(has_reads = T)
for(i in 1:nrow(seeds.tbl)){
  cat('  collecting reads for', seeds.tbl\$Header[i], '...')
  ss_files <- list.files(tmp_folder, pattern = paste0('_', seeds.tbl\$Header[i], '_'), full.names = T)
  cat('found', length(ss_files), 'such files\n')
  if(length(ss_files) > 0){
    ss1_files <- sort(ss_files[grep('R1', ss_files)])
    seed1_file <- file.path(tmp_folder, paste0(seeds.tbl\$Header[i], '_R1.fastq'))
    if(file.exists(seed1_file)) ok <- file.remove(seed1_file)
    ok <- file.append(seed1_file, ss1_files)
    ss2_files <- sort(ss_files[grep('R2', ss_files)])
    seed2_file <- file.path(tmp_folder, paste0(seeds.tbl\$Header[i], '_R2.fastq'))
    if(file.exists(seed2_file)) ok <- file.remove(seed2_file)
    ok <- file.append(seed2_file, ss2_files)
    ok <- writeFasta(seeds.tbl[i,], out.file = file.path(tmp_folder, paste0(seeds.tbl\$Header[i], '.fasta')))
  } else {
    seeds.tbl\$has_reads[i] <- F
  }
}
seeds.tbl <- filter(seeds.tbl, has_reads)
cat('Found reads for', nrow(seeds.tbl), 'seeds\n')
writeLines(seeds.tbl\$Header, con = file.path(tmp_folder, 'seed_names.txt'))" $seed_sequence_file $tmp_folder


echo "*** The seed assemblies:"
assembly_folder=$tmp_folder/assembly
if [ -d $assembly_folder ]
then
  rm -rf $assembly_folder
fi
seed_names=$(cat $tmp_folder/seed_names.txt)
for seed in $seed_names
do
  echo "*** Assembly for seed $seed"
  if [ ! -d $assembly_folder ]
  then
    mkdir -p $assembly_folder
  fi
  r1=$tmp_folder/$seed\_R1.fastq
  r2=$tmp_folder/$seed\_R2.fastq
  fa=$tmp_folder/$seed.fasta
  $spades_app spades.py -o $assembly_folder --tmp-dir $assembly_folder --threads $threads \
  --isolate --trusted-contigs $fa --pe-1 1 $r1 --pe-2 1 $r2 > $assembly_folder/screen.log
  if [ -f $assembly_folder/contigs.fasta ]
  then
    $blast_app blastn -db $dbase -query $assembly_folder/contigs.fasta -num_threads $threads \
      -outfmt '6 qseqid sseqid pident qlen qstart qend slen sstart send bitscore' \
      -penalty -1 -reward 1 -gapopen 2 -gapextend 1 -out $out_folder/blst.txt
    Rscript -e "suppressMessages(library(tidyverse))
      suppressMessages(library(microseq))
      args <- commandArgs(trailingOnly = T)
      out_folder <- args[1]
      assembly_folder <- args[2]
      seed_id <- args[3]
      seed.tbl <- readFasta(file.path(assembly_folder, 'contigs.fasta')) %>% 
        mutate(length = str_length(Sequence)) %>% 
        slice_max(length)
      b.tbl <- suppressMessages(read_delim(file.path(out_folder, 'blst.txt'), delim = '\t',
                          col_names = c('qseqid', 'sseqid', 'pident', 'qlen',
                                        'qstart', 'qend', 'slen', 'sstart', 'send', 'bitscore'))) %>% 
        filter(qseqid %in% seed.tbl\$Header) %>% 
        arrange(desc(bitscore)) %>% 
        slice(1)
      if(b.tbl\$sstart > b.tbl\$send){
        seed.tbl <- seed.tbl %>% mutate(Sequence = reverseComplement(Sequence))
      }
      seed.tbl <- seed.tbl %>%
        mutate(Header = str_c(seed_id, ' best_match=', b.tbl\$sseqid[1], ';identity=', b.tbl\$pident[1]))
      writeFasta(seed.tbl, out.file = file.path(out_folder, str_c(seed_id, '.fasta')))
      " $out_folder $assembly_folder $seed
  fi
  rm -rf $assembly_folder
done
rm $out_folder/blst.txt

