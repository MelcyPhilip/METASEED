#!/bin/bash
seed_file=$1
seed_reads_folder=$2
tmp_folder=$3
contigs_folder=$4
threads=$5


##############
### Settings
###
memory=50
assembly_folder=$tmp_folder/assembly
if [ -d $assembly_folder ]
then
  rm -rf $assembly_folder
fi


#####################################
### Collecting reads for each seed
###
echo "**********************************************"
echo "Collecting reads for each seed..."
module load R/4.1.0
Rscript -e "args <- commandArgs(trailingOnly = T)
library(microseq)
seed_file <- args[1]
seed_reads_folder <- args[2]
tmp_folder <- args[3]

seeds.tbl <- readFasta(seed_file) %>% 
  mutate(Header = word(Header, 1))
for(i in 1:nrow(seeds.tbl)){
  cat('  collecting reads for', seeds.tbl\$Header[i], '...\n')
  s_files <- list.files(seed_reads_folder, pattern = paste0('_', seeds.tbl\$Header[i], '_'), full.names = T)
  if(length(s_files) > 0){
    s1_files <- sort(s_files[grep('R1', s_files)])
    s_file <- file.path(seed_reads_folder, paste0(seeds.tbl\$Header[i], '_R1.fastq'))
    if(file.exists(s_file)) ok <- file.remove(s_file)
    ok <- file.append(s_file, s1_files)
    s2_files <- sort(s_files[grep('R2', s_files)])
    s_file <- file.path(seed_reads_folder, paste0(seeds.tbl\$Header[i], '_R2.fastq'))
    if(file.exists(s_file)) ok <- file.remove(s_file)
    ok <- file.append(s_file, s2_files)
    ok <- writeFasta(seeds.tbl[i,], out.file = file.path(seed_reads_folder, paste0(seeds.tbl\$Header[i], '.fasta')))
  }
}
writeLines(seeds.tbl\$Header, con = file.path(tmp_folder, 'seed_names.txt'))
cat('done!\n')" $seed_file $seed_reads_folder $tmp_folder
module purge


############################
### The assemblies
###
echo "*****************************************"
echo "The seed assemblies"
seed_names=$(cat $tmp_folder/seed_names.txt)
for seed in $seed_names
do
  if [ ! -d $assembly_folder ]
  then
    mkdir -p $assembly_folder
  fi
  r1=$seed_reads_folder/$seed\_R1.fastq
  r2=$seed_reads_folder/$seed\_R2.fastq
  fa=$seed_reads_folder/$seed.fasta
  echo "***********************************************"
  echo "Assembly of $seed"

  singularity exec /cvmfs/singularity.galaxyproject.org/s/p/spades:3.15.3--h95f258a_0 spades.py \
  -o $assembly_folder \
  --tmp-dir $assembly_folder \
  --memory $memory \
  --threads $threads \
  --isolate \
  --trusted-contigs $fa \
  --pe-1 1 $r1 \
  --pe-2 1 $r2 \
  > $assembly_folder/screen.log
  
  if [ -f $assembly_folder/contigs.fasta ]
  then
    echo "  successful!"
    cp $assembly_folder/contigs.fasta $contigs_folder/$seed\_contigs.fasta
  fi
  rm -rf $assembly_folder
done
echo "********************************"
echo "done!"
