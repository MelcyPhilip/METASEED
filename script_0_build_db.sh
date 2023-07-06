#!/bin/bash
seed_file=$1
dbase=$2


################
### Settings
###
blast_exe="singularity exec /METASEED/blast:2.12.0--hf3cf87c_4.sif" #path to blast executable
bowtie2_exe="singularity exec /METASEED/bowtie2:2.5.0--py39h3321a2d_0.sif" #path to bowtie2 executable


###########################################################
### Building the bowtie index
###
echo "*************************************************************"
echo "Building bowtie2 index from $seed_file"
$bowtie2_exe bowtie2-build \
 $seed_file \
 $dbase


###########################################################
### Building the blast database
###
echo "*************************************************************"
echo "Building blast database from $seed_file"
$blast_exe makeblastdb \
 -dbtype nucl \
 -in $seed_file \
 -out $dbase
