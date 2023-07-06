#!/bin/bash

#SBATCH --array=1-24                     #Use your prefered number of array jobs
#SBATCH --nodes=1                        # We always use 1 node
#SBATCH --ntasks=10                      # The number of threads reserved
#SBATCH --mem=50G                        # The amount of memory reserved
#SBATCH --partition=smallmem,hugemem     # For < 100GB use smallmem, for >100GB use hugemem
#SBATCH --time=24:00:00                  # Runs for maximum this time
#SBATCH --job-name=seed24_assembly         # Sensible name for the job
#SBATCH --output=seed24_assembly_%j_%a.log # Logfile output here
#SBATCH --exclude=cn-14

###############
### Settings
##
seed_file=/path to the seed file/seedfile.fasta
seed_abundance_file=/path to the seed readcount file/seedreadcount.txt
fastq_folder=/path to shotgun metagenome raw data folder/shotgundata_folder
metadata_file=/path to metadata file/metadata_file.txt
work_folder=work
tmp_folder=$work_folder/tmp
aln_folder=$work_folder/aln
seed_reads_folder=$work_folder/seed_reads
contigs_folder=$work_folder/rec_16S
dbase=$tmp_folder/seeds
threads=10 # The number of threads reserved

##########################
### Initialization
###
if [ ! -d $work_folder ]
then
  mkdir -p $work_folder
fi
if [ ! -d $tmp_folder ]
then
  mkdir -p $tmp_folder
fi
if [ ! -d $aln_folder ]
then
  mkdir -p $aln_folder
fi
if [ ! -d $seed_reads_folder ]
then
  mkdir -p $seed_reads_folder
fi
if [ ! -d $contigs_folder ]
then
  mkdir -p $contigs_folder
fi
SampleID=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++){f[$i] = i}}{print $(f["SampleID"])}' $metadata_file)
Rawfile_R1=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++){f[$i] = i}}{print $(f["Rawfile_R1"])}' $metadata_file)
Rawfile_R2=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++){f[$i] = i}}{print $(f["Rawfile_R2"])}' $metadata_file)
line=$(($SLURM_ARRAY_TASK_ID+1))
prefix=$(echo $SampleID | awk -vidx=$line '{print $idx}')
raw_r1=$fastq_folder/$(echo $Rawfile_R1 | awk -vidx=$line '{print $idx}')
raw_r2=$fastq_folder/$(echo $Rawfile_R2 | awk -vidx=$line '{print $idx}')

##################################################
### Building BLAST database from seed sequences
### NB! No array jobs here!
###
# ./script_0_build_db.sh \
#  $seed_file \
#  $dbase


#################################################
### First bowtie2 mapping, then blast alignment
###
# ./script_1_align.sh \
#  $prefix \
#  $raw_r1 \
#  $raw_r2 \
#  $dbase \
#  $aln_folder \
#  $threads
  

##################################
### Collecting the seed reads
## 
 # module load R/4.1.0
 # Rscript script_2_collect_reads.R \
 # $prefix \
 # $seed_abundance_file \
 # $raw_r1 \
 # $raw_r2 \
 # $aln_folder \
 # $seed_reads_folder \

##################################
### The assembly
### NB! No array jobs here!
###
#./script_3_assembly.sh \
# $seed_file \
# $seed_reads_folder \
# $tmp_folder \
# $contigs_folder \
# $threads


#############################
### Collecting estimated 16S
### NB! No array jobs here!
###
   # module load R/4.1.0
   # Rscript script_4_finalize.R \
   # $contigs_folder


