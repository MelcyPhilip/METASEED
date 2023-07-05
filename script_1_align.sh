#!/bin/bash
prefix=$1
r1=$2
r2=$3
dbase=$4
aln_folder=$5
threads=$6


####################
### Settings
###
blast_exe="singularity exec /METASEED/blast:2.12.0--hf3cf87c_4.sif"
bowtie2_exe="singularity exec /METASEED/bowtie2:2.5.0--py39h3321a2d_0.sif"
samtools_exe=="singularity exec /METASEED/samtools:1.9--h91753b0_8.sif"


###############
### Out files
###
r1_fasta=$aln_folder/$prefix\_R1.fasta
r2_fasta=$aln_folder/$prefix\_R2.fasta
r1_txt=$aln_folder/$prefix\_R1.txt
r2_txt=$aln_folder/$prefix\_R2.txt


###################
### Readmapping
###
echo "**************************************************************************"
echo "Mapping $r1"
echo "Writing to $r1_fasta"
$bowtie2_exe bowtie2 \
 --threads $threads \
 --local \
 -x $dbase \
 -U $r1 \
 | $samtools_exe samtools view \
 --threads $threads -b - \
 | $samtools_exe samtools fasta \
 -0 $r1_fasta -F 4 -

echo "*************************************************************************"
echo "Mapping $r2"
echo "Writing to $r2_fasta"
$bowtie2_exe bowtie2 \
 --threads $threads \
 --local \
 -x $dbase \
 -U $r2 \
 | $samtools_exe samtools view \
 --threads $threads -b - \
 | $samtools_exe samtools fasta \
 -0 $r2_fasta -F 4 -
 
 
 
#######################
### Aligning
###
echo "*************************************************************************"
echo "Aligning $r1_fasta"
echo "Writing to $r1_txt"
blast_exe \
  -db $dbase \
  -query $r1_fasta \
  -num_alignments 250 \
  -perc_identity 99 \
  -num_threads $threads \
  -outfmt '6 qseqid sseqid pident gapopen qlen qstart qend slen sstart send bitscore' \
  -penalty -1 -reward 1 -gapopen 2 -gapextend 1 \
  -out $r1_txt

echo "*************************************************************************"
echo "Aligning $r2_fasta"
echo "Writing to $r2_txt"
blast_exe \
  -db $dbase \
  -query $r2_fasta \
  -num_alignments 250 \
  -perc_identity 99 \
  -num_threads $threads \
  -outfmt '6 qseqid sseqid pident gapopen qlen qstart qend slen sstart send bitscore' \
  -penalty -1 -reward 1 -gapopen 2 -gapextend 1 \
  -out $r2_txt


###################
### Cleaning
###
rm $r1_fasta
rm $r2_fasta

