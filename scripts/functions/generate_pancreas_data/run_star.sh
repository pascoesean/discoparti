#!/usr/bin/env bash

#SBATCH -c 4		# number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#SBATCH --exclude=c17,c21 # this node is ass
#SBATCH --exclusive
#SBATCH --job-name STAR_alignment-HPAP140-160	# Job name
#SBATCH --output=R-%x.%j.out		# File to which standard out will be written
#SBATCH --error=R-%x.%j.err 		# File to which standard err will be written  # /home/Genomes/STAR_indexes/hg38_star_topLevel
#SBATCH --mail-type=ALL
#SBATCH --mail-user=spascoe@mit.edu

DONORS=("055 064 123 130 135 137 138 141 156 167 172 131 155 082 140 024 029 037 049 146 036 080")

module load star/2.7.9a

STAR --genomeLoad LoadAndExit --genomeDir reference_genome/hg38_index

for donor in $DONORS; do

    R1_FILES=($(ls hpapdata/HPAP-$(basename $donor)/fastq/*R1*fastq.gz))
    R2_FILES=($(ls hpapdata/HPAP-$(basename $donor)/fastq/*R2*fastq.gz))


    for i in ${!R1_FILES[@]}; do 

        echo "${R1_FILES[i]} + ${R2_FILES[i]}"

        STAR --genomeDir reference_genome/hg38_index \
             --genomeLoad LoadAndKeep \
             --runThreadN 4 \
             --readFilesCommand zcat \
             --readFilesIn ${R2_FILES[i]} ${R1_FILES[i]} \
             --outFileNamePrefix hpapdata/HPAP-$(basename $donor)/STAR-${i} \
             --soloType CB_UMI_Simple \
             --soloCBwhitelist reference_genome/x10_3p3v_whitelist.txt \
             --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
             --soloCellFilter EmptyDrops_CR \
             --soloOutFileNames  Solo.out/ genes.tsv barcodes.tsv matrix.mtx \
             --limitBAMsortRAM 10000000000 \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMunmapped Within \
             --outSAMattributes Standard
    done

    echo " **FINISHED DONOR: $(basename $donor) **. whoopie!"

done

STAR --genomeLoad Remove --genomeDir reference_genome/hg38_index