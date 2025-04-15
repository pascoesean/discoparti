#!/usr/bin/env bash

#SBATCH -c 6		# number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#SBATCH --exclude=c17 # this node is ass
#SBATCH --exclusive
#SBATCH --job-name generate-fastqc-064-123		# Job name
#SBATCH -o %j.out			# File to which standard out will be written
#SBATCH -e %j.err 		# File to which standard err will be written  # /home/Genomes/STAR_indexes/hg38_star_topLevel
#SBATCH --mail-type=ALL
#SBATCH --mail-user=spascoe@mit.edu

DONORS=("064 072 074 075 077 080 082 092 093 095 097 099 101 103 105 110 114 117 118 119 122 123")


module load fastqc

for donor in $DONORS; do

     mkdir hpapdata/HPAP-$(basename $donor)/fastqc_results

     R2_FILES=$(ls hpapdata/HPAP-$(basename $donor)/fastq/*R2*fastq.gz)

     echo "IM DOING: $R2_FILES"

     fastqc $R2_FILES --outdir hpapdata/HPAP-$(basename $donor)/fastqc_results --threads 6

    
    echo " **FINISHED DONOR: $(basename $donor) **. whoopie!"

done
