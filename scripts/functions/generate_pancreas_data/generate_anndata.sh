#!/usr/bin/env bash

#SBATCH --exclude=c17 # this node is ass
#SBATCH --exclusive
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --job-name generate-anndataobj		# Job name
#SBATCH -o %j.out			# File to which standard out will be written
#SBATCH -e %j.err 		# File to which standard err will be written  # /home/Genomes/STAR_indexes/hg38_star_topLevel
#SBATCH --mail-type=ALL
#SBATCH --mail-user=spascoe@mit.edu

source "/net/bmc-lab8/data/lab/dallgglab/users/spascoe/conda/etc/profile.d/conda.sh" 
source "/net/bmc-lab8/data/lab/dallgglab/users/spascoe/conda/etc/profile.d/mamba.sh"
mamba activate pancre-yas-env

python3 discoparti/scripts/functions/generate_pancreas_data/generate_anndata_obj.py