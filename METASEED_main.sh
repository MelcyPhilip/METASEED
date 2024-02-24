#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --job-name=metaseed
#SBATCH --output=metaseed_%j.log


###############
### Settings
##
seed_sequence_file="example/seed_sequences.fasta"   # fasta file with the seeds (=amplicon sequences)
seed_abundance_file="example/seed_abundances.tsv"   # text file with abundances for each seed
sample_table_file="example/sample_data.tsv"         # text file with metadat for each sample
tmp_folder="tmp"                                    # folder for temporary output
out_folder="out"                                    # folder for final results
threads=10                                          # number of threads
dbase=$tmp_folder/METASEED                          # the bowtie2 and BLASt database prefix
max_n=3                                             # maximum number of seeds a read can match


###############################
### The software containers
###
export bowtie2_app="apptainer exec containers/bowtie2:2.5.0--py39h3321a2d_0.sif"
export blast_app="apptainer exec containers/blast:2.12.0--hf3cf87c_4.sif"
export samtools_app="apptainer exec containers/samtools:1.9--h91753b0_8.sif"
export bbmap_app="apptainer exec containers/bbmap%3A39.06--h92535d8_0.sif"
export spades_app="apptainer exec containers/spades:3.15.3--h95f258a_0.sif"


#######################
### Initializing
###
if [ ! -d $tmp_folder ]
then
  mkdir -p $tmp_folder
fi
if [ ! -d $out_folder ]
then
  mkdir -p $out_folder
fi
module load R/4.3.1     # Make certain R is available. This is how we load R on our HPC cluster


echo "************************************************************"
echo "************************************************************"
echo "*** Step 1"
echo "*** Creates the databases for the alignments"
$bowtie2_app bowtie2-build $seed_sequence_file $dbase
$blast_app makeblastdb -dbtype nucl -in $seed_sequence_file -out $dbase


echo
echo
echo "************************************************************"
echo "************************************************************"
echo "*** Step 2"
echo "*** Collecting reads from each sample. This takes time!"
echo "*** Looping over samples..."
sample_id_column=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++){f[$i] = i}}{print $(f["sample_id"])}' $sample_table_file)
R1_file_column=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++){f[$i] = i}}{print $(f["R1_file"])}' $sample_table_file)
R2_file_column=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++){f[$i] = i}}{print $(f["R2_file"])}' $sample_table_file)
n_rows=$(wc $sample_table_file | awk -vj=1 '{print $j}')
for i in $(seq 2 $n_rows)
do
  sample_id=$(echo $sample_id_column | awk -vj=$i '{print $j}')
  R1_file=$(echo $R1_file_column | awk -vj=$i '{print $j}')
  R2_file=$(echo $R2_file_column | awk -vj=$i '{print $j}')
  echo "************************************************************"
  echo "*** Sample $sample_id"
  ./read_collecting.sh $sample_id $R1_file $R2_file $dbase $tmp_folder $seed_abundance_file $max_n $threads
done


echo
echo
echo "************************************************************"
echo "************************************************************"
echo "*** Step 3"
echo "*** Assembly of reads for each seed"
./assembly.sh $seed_sequence_file $tmp_folder $out_folder $threads $dbase




