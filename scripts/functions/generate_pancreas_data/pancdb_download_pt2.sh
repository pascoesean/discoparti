#!/usr/bin/env bash
#SBATCH --job-name=get-pancreas-data-pt2
#SBATCH --exclude=c17,c21 # this node is ass
#SBATCH -N 1                    ## Number of Nodes
#SBATCH --mail-type=ALL
#SBATCH --mail-user=spascoe@mit.edu


#
# This script downloads HPAP data/metadata into your working directory:
#   * Data (if requested) will be saved in "hpapdata" subdirectory;
#   * Metadata (if requested) will be saved in "metadata" subdirectory.

DATA_SERVER="https://hpapdata.faryabilab.com"

FILES="
HPAP-128/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-128_GEX_FGC2493_112654_S29_L001_I1_001.fastq.gz
HPAP-128/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-128_GEX_FGC2493_112654_S29_L001_I2_001.fastq.gz
HPAP-128/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-128_GEX_FGC2493_112654_S29_L001_R1_001.fastq.gz
HPAP-128/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-128_GEX_FGC2493_112654_S29_L001_R2_001.fastq.gz
HPAP-129/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-129_GEX_FGC2493_112656_S7_L002_I1_001.fastq.gz
HPAP-129/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-129_GEX_FGC2493_112656_S7_L002_I2_001.fastq.gz
HPAP-129/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-129_GEX_FGC2493_112656_S7_L002_R1_001.fastq.gz
HPAP-129/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-129_GEX_FGC2493_112656_S7_L002_R2_001.fastq.gz
HPAP-129/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-129_GEX_FGC2578_112656_S12_L001_I1_001.fastq.gz
HPAP-129/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-129_GEX_FGC2578_112656_S12_L001_I2_001.fastq.gz
HPAP-129/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-129_GEX_FGC2578_112656_S12_L001_R1_001.fastq.gz
HPAP-129/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-129_GEX_FGC2578_112656_S12_L001_R2_001.fastq.gz
HPAP-129/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-129_GEX_FGC2578_112656_S12_L002_I1_001.fastq.gz
HPAP-129/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-129_GEX_FGC2578_112656_S12_L002_I2_001.fastq.gz
HPAP-129/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-129_GEX_FGC2578_112656_S12_L002_R1_001.fastq.gz
HPAP-129/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-129_GEX_FGC2578_112656_S12_L002_R2_001.fastq.gz
HPAP-130/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-130_GEX_FGC2493_112658_S30_L001_I1_001.fastq.gz
HPAP-130/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-130_GEX_FGC2493_112658_S30_L001_I2_001.fastq.gz
HPAP-130/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-130_GEX_FGC2493_112658_S30_L001_R1_001.fastq.gz
HPAP-130/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-130_GEX_FGC2493_112658_S30_L001_R2_001.fastq.gz
HPAP-130/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-130_GEX_FGC2578_112658_S13_L001_I1_001.fastq.gz
HPAP-130/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-130_GEX_FGC2578_112658_S13_L001_I2_001.fastq.gz
HPAP-130/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-130_GEX_FGC2578_112658_S13_L001_R1_001.fastq.gz
HPAP-130/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-130_GEX_FGC2578_112658_S13_L001_R2_001.fastq.gz
HPAP-130/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-130_GEX_FGC2578_112658_S13_L002_I1_001.fastq.gz
HPAP-130/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-130_GEX_FGC2578_112658_S13_L002_I2_001.fastq.gz
HPAP-130/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-130_GEX_FGC2578_112658_S13_L002_R1_001.fastq.gz
HPAP-130/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-130_GEX_FGC2578_112658_S13_L002_R2_001.fastq.gz
HPAP-131/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-131_GEX_FGC2493_112660_S31_L001_I1_001.fastq.gz
HPAP-131/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-131_GEX_FGC2493_112660_S31_L001_I2_001.fastq.gz
HPAP-131/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-131_GEX_FGC2493_112660_S31_L001_R1_001.fastq.gz
HPAP-131/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-131_GEX_FGC2493_112660_S31_L001_R2_001.fastq.gz
HPAP-131/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-131_GEX_FGC2578_112660_S14_L001_I1_001.fastq.gz
HPAP-131/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-131_GEX_FGC2578_112660_S14_L001_I2_001.fastq.gz
HPAP-131/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-131_GEX_FGC2578_112660_S14_L001_R1_001.fastq.gz
HPAP-131/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-131_GEX_FGC2578_112660_S14_L001_R2_001.fastq.gz
HPAP-131/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-131_GEX_FGC2578_112660_S14_L002_I1_001.fastq.gz
HPAP-131/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-131_GEX_FGC2578_112660_S14_L002_I2_001.fastq.gz
HPAP-131/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-131_GEX_FGC2578_112660_S14_L002_R1_001.fastq.gz
HPAP-131/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-131_GEX_FGC2578_112660_S14_L002_R2_001.fastq.gz
HPAP-135/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-135_GEX_FGC2578_113766_S2_L001_I1_001.fastq.gz
HPAP-135/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-135_GEX_FGC2578_113766_S2_L001_I2_001.fastq.gz
HPAP-135/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-135_GEX_FGC2578_113766_S2_L001_R1_001.fastq.gz
HPAP-135/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-135_GEX_FGC2578_113766_S2_L001_R2_001.fastq.gz
HPAP-135/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-135_GEX_FGC2578_113766_S2_L002_I1_001.fastq.gz
HPAP-135/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-135_GEX_FGC2578_113766_S2_L002_I2_001.fastq.gz
HPAP-135/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-135_GEX_FGC2578_113766_S2_L002_R1_001.fastq.gz
HPAP-135/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-135_GEX_FGC2578_113766_S2_L002_R2_001.fastq.gz
HPAP-136/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-136_GEX_FGC2578_113767_S3_L001_I1_001.fastq.gz
HPAP-136/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-136_GEX_FGC2578_113767_S3_L001_I2_001.fastq.gz
HPAP-136/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-136_GEX_FGC2578_113767_S3_L001_R1_001.fastq.gz
HPAP-136/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-136_GEX_FGC2578_113767_S3_L001_R2_001.fastq.gz
HPAP-136/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-136_GEX_FGC2578_113767_S3_L002_I1_001.fastq.gz
HPAP-136/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-136_GEX_FGC2578_113767_S3_L002_I2_001.fastq.gz
HPAP-136/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-136_GEX_FGC2578_113767_S3_L002_R1_001.fastq.gz
HPAP-136/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-136_GEX_FGC2578_113767_S3_L002_R2_001.fastq.gz
HPAP-137/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-137_GEX_FGC2578_113768_S4_L001_I1_001.fastq.gz
HPAP-137/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-137_GEX_FGC2578_113768_S4_L001_I2_001.fastq.gz
HPAP-137/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-137_GEX_FGC2578_113768_S4_L001_R1_001.fastq.gz
HPAP-137/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-137_GEX_FGC2578_113768_S4_L001_R2_001.fastq.gz
HPAP-137/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-137_GEX_FGC2578_113768_S4_L002_I1_001.fastq.gz
HPAP-137/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-137_GEX_FGC2578_113768_S4_L002_I2_001.fastq.gz
HPAP-137/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-137_GEX_FGC2578_113768_S4_L002_R1_001.fastq.gz
HPAP-137/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-137_GEX_FGC2578_113768_S4_L002_R2_001.fastq.gz
HPAP-138/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-138_GEX_FGC2578_113769_S5_L001_I1_001.fastq.gz
HPAP-138/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-138_GEX_FGC2578_113769_S5_L001_I2_001.fastq.gz
HPAP-138/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-138_GEX_FGC2578_113769_S5_L001_R1_001.fastq.gz
HPAP-138/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-138_GEX_FGC2578_113769_S5_L001_R2_001.fastq.gz
HPAP-138/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-138_GEX_FGC2578_113769_S5_L002_I1_001.fastq.gz
HPAP-138/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-138_GEX_FGC2578_113769_S5_L002_I2_001.fastq.gz
HPAP-138/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-138_GEX_FGC2578_113769_S5_L002_R1_001.fastq.gz
HPAP-138/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-138_GEX_FGC2578_113769_S5_L002_R2_001.fastq.gz
HPAP-139/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-139_GEX_FGC2578_113770_S6_L001_I1_001.fastq.gz
HPAP-139/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-139_GEX_FGC2578_113770_S6_L001_I2_001.fastq.gz
HPAP-139/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-139_GEX_FGC2578_113770_S6_L001_R1_001.fastq.gz
HPAP-139/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-139_GEX_FGC2578_113770_S6_L001_R2_001.fastq.gz
HPAP-139/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-139_GEX_FGC2578_113770_S6_L002_I1_001.fastq.gz
HPAP-139/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-139_GEX_FGC2578_113770_S6_L002_I2_001.fastq.gz
HPAP-139/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-139_GEX_FGC2578_113770_S6_L002_R1_001.fastq.gz
HPAP-139/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-139_GEX_FGC2578_113770_S6_L002_R2_001.fastq.gz
HPAP-140/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-140_GEX_FGC2578_113774_S8_L001_I1_001.fastq.gz
HPAP-140/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-140_GEX_FGC2578_113774_S8_L001_I2_001.fastq.gz
HPAP-140/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-140_GEX_FGC2578_113774_S8_L001_R1_001.fastq.gz
HPAP-140/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-140_GEX_FGC2578_113774_S8_L001_R2_001.fastq.gz
HPAP-140/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-140_GEX_FGC2578_113774_S8_L002_I1_001.fastq.gz
HPAP-140/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-140_GEX_FGC2578_113774_S8_L002_I2_001.fastq.gz
HPAP-140/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-140_GEX_FGC2578_113774_S8_L002_R1_001.fastq.gz
HPAP-140/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-140_GEX_FGC2578_113774_S8_L002_R2_001.fastq.gz
HPAP-141/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-141_GEX_FGC2578_113771_S7_L001_I1_001.fastq.gz
HPAP-141/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-141_GEX_FGC2578_113771_S7_L001_I2_001.fastq.gz
HPAP-141/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-141_GEX_FGC2578_113771_S7_L001_R1_001.fastq.gz
HPAP-141/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-141_GEX_FGC2578_113771_S7_L001_R2_001.fastq.gz
HPAP-141/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-141_GEX_FGC2578_113771_S7_L002_I1_001.fastq.gz
HPAP-141/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-141_GEX_FGC2578_113771_S7_L002_I2_001.fastq.gz
HPAP-141/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-141_GEX_FGC2578_113771_S7_L002_R1_001.fastq.gz
HPAP-141/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-141_GEX_FGC2578_113771_S7_L002_R2_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_FGC2566_113871_S3_L001_I1_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_FGC2566_113871_S3_L001_I2_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_FGC2566_113871_S3_L001_R1_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_FGC2566_113871_S3_L001_R2_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_FGC2566_113871_S3_L002_I1_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_FGC2566_113871_S3_L002_I2_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_FGC2566_113871_S3_L002_R1_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_FGC2566_113871_S3_L002_R2_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_FGC2569_113871_S1_L001_I1_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_FGC2569_113871_S1_L001_I2_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_FGC2569_113871_S1_L001_R1_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_FGC2569_113871_S1_L001_R2_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_FGC2569_113871_S1_L002_I1_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_FGC2569_113871_S1_L002_I2_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_FGC2569_113871_S1_L002_R1_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_FGC2569_113871_S1_L002_R2_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_NSS0018_113871_S13_L001_I1_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_NSS0018_113871_S13_L001_I2_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_NSS0018_113871_S13_L001_R1_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_NSS0018_113871_S13_L001_R2_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_NSS0018_113871_S13_L002_I1_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_NSS0018_113871_S13_L002_I2_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_NSS0018_113871_S13_L002_R1_001.fastq.gz
HPAP-146/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-146_GEX_NSS0018_113871_S13_L002_R2_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0018_5000117_S6_L001_I1_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0018_5000117_S6_L001_I2_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0018_5000117_S6_L001_R1_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0018_5000117_S6_L001_R2_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0018_5000117_S6_L002_I1_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0018_5000117_S6_L002_I2_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0018_5000117_S6_L002_R1_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0018_5000117_S6_L002_R2_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0020_5000117_S10_L001_I1_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0020_5000117_S10_L001_I2_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0020_5000117_S10_L001_R1_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0020_5000117_S10_L001_R2_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0020_5000117_S10_L002_I1_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0020_5000117_S10_L002_I2_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0020_5000117_S10_L002_R1_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0020_5000117_S10_L002_R2_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0025_5000117_S12_L001_I1_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0025_5000117_S12_L001_I2_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0025_5000117_S12_L001_R1_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0025_5000117_S12_L001_R2_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0025_5000117_S12_L002_I1_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0025_5000117_S12_L002_I2_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0025_5000117_S12_L002_R1_001.fastq.gz
HPAP-155/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-155_GEX_NSS0025_5000117_S12_L002_R2_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0018_5000118_S7_L001_I1_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0018_5000118_S7_L001_I2_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0018_5000118_S7_L001_R1_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0018_5000118_S7_L001_R2_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0018_5000118_S7_L002_I1_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0018_5000118_S7_L002_I2_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0018_5000118_S7_L002_R1_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0018_5000118_S7_L002_R2_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0020_5000118_S11_L001_I1_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0020_5000118_S11_L001_I2_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0020_5000118_S11_L001_R1_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0020_5000118_S11_L001_R2_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0020_5000118_S11_L002_I1_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0020_5000118_S11_L002_I2_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0020_5000118_S11_L002_R1_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0020_5000118_S11_L002_R2_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0025_5000118_S8_L001_I1_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0025_5000118_S8_L001_I2_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0025_5000118_S8_L001_R1_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0025_5000118_S8_L001_R2_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0025_5000118_S8_L002_I1_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0025_5000118_S8_L002_I2_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0025_5000118_S8_L002_R1_001.fastq.gz
HPAP-156/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-156_GEX_NSS0025_5000118_S8_L002_R2_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0018_5000119_S8_L001_I1_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0018_5000119_S8_L001_I2_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0018_5000119_S8_L001_R1_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0018_5000119_S8_L001_R2_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0018_5000119_S8_L002_I1_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0018_5000119_S8_L002_I2_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0018_5000119_S8_L002_R1_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0018_5000119_S8_L002_R2_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0020_5000119_S12_L001_I1_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0020_5000119_S12_L001_I2_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0020_5000119_S12_L001_R1_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0020_5000119_S12_L001_R2_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0020_5000119_S12_L002_I1_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0020_5000119_S12_L002_I2_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0020_5000119_S12_L002_R1_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0020_5000119_S12_L002_R2_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0025_5000119_S11_L001_I1_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0025_5000119_S11_L001_I2_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0025_5000119_S11_L001_R1_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0025_5000119_S11_L001_R2_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0025_5000119_S11_L002_I1_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0025_5000119_S11_L002_I2_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0025_5000119_S11_L002_R1_001.fastq.gz
HPAP-157/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-157_GEX_NSS0025_5000119_S11_L002_R2_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0018_5000122_S11_L001_I1_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0018_5000122_S11_L001_I2_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0018_5000122_S11_L001_R1_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0018_5000122_S11_L001_R2_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0018_5000122_S11_L002_I1_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0018_5000122_S11_L002_I2_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0018_5000122_S11_L002_R1_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0018_5000122_S11_L002_R2_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0020_5000122_S13_L001_I1_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0020_5000122_S13_L001_I2_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0020_5000122_S13_L001_R1_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0020_5000122_S13_L001_R2_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0020_5000122_S13_L002_I1_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0020_5000122_S13_L002_I2_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0020_5000122_S13_L002_R1_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0020_5000122_S13_L002_R2_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0025_5000122_S3_L001_I1_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0025_5000122_S3_L001_I2_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0025_5000122_S3_L001_R1_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0025_5000122_S3_L001_R2_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0025_5000122_S3_L002_I1_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0025_5000122_S3_L002_I2_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0025_5000122_S3_L002_R1_001.fastq.gz
HPAP-160/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-160_GEX_NSS0025_5000122_S3_L002_R2_001.fastq.gz
HPAP-167/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-167_GEX_NSS0020_5000237_S5_L001_I1_001.fastq.gz
HPAP-167/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-167_GEX_NSS0020_5000237_S5_L001_I2_001.fastq.gz
HPAP-167/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-167_GEX_NSS0020_5000237_S5_L001_R1_001.fastq.gz
HPAP-167/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-167_GEX_NSS0020_5000237_S5_L001_R2_001.fastq.gz
HPAP-167/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-167_GEX_NSS0020_5000237_S5_L002_I1_001.fastq.gz
HPAP-167/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-167_GEX_NSS0020_5000237_S5_L002_I2_001.fastq.gz
HPAP-167/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-167_GEX_NSS0020_5000237_S5_L002_R1_001.fastq.gz
HPAP-167/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-167_GEX_NSS0020_5000237_S5_L002_R2_001.fastq.gz
HPAP-167/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-167_GEX_NSS0025_5000237_S7_L001_I1_001.fastq.gz
HPAP-167/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-167_GEX_NSS0025_5000237_S7_L001_I2_001.fastq.gz
HPAP-167/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-167_GEX_NSS0025_5000237_S7_L001_R1_001.fastq.gz
HPAP-167/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-167_GEX_NSS0025_5000237_S7_L001_R2_001.fastq.gz
HPAP-167/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-167_GEX_NSS0025_5000237_S7_L002_I1_001.fastq.gz
HPAP-167/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-167_GEX_NSS0025_5000237_S7_L002_I2_001.fastq.gz
HPAP-167/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-167_GEX_NSS0025_5000237_S7_L002_R1_001.fastq.gz
HPAP-167/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-167_GEX_NSS0025_5000237_S7_L002_R2_001.fastq.gz
HPAP-169/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-169_GEX_NSS0020_5000239_S7_L001_I1_001.fastq.gz
HPAP-169/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-169_GEX_NSS0020_5000239_S7_L001_I2_001.fastq.gz
HPAP-169/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-169_GEX_NSS0020_5000239_S7_L001_R1_001.fastq.gz
HPAP-169/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-169_GEX_NSS0020_5000239_S7_L001_R2_001.fastq.gz
HPAP-169/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-169_GEX_NSS0020_5000239_S7_L002_I1_001.fastq.gz
HPAP-169/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-169_GEX_NSS0020_5000239_S7_L002_I2_001.fastq.gz
HPAP-169/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-169_GEX_NSS0020_5000239_S7_L002_R1_001.fastq.gz
HPAP-169/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-169_GEX_NSS0020_5000239_S7_L002_R2_001.fastq.gz
HPAP-171/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-171_GEX_NSS0032_5000904_S1_L001_I1_001.fastq.gz
HPAP-171/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-171_GEX_NSS0032_5000904_S1_L001_I2_001.fastq.gz
HPAP-171/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-171_GEX_NSS0032_5000904_S1_L001_R1_001.fastq.gz
HPAP-171/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-171_GEX_NSS0032_5000904_S1_L001_R2_001.fastq.gz
HPAP-171/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-171_GEX_NSS0032_5000904_S1_L002_I1_001.fastq.gz
HPAP-171/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-171_GEX_NSS0032_5000904_S1_L002_I2_001.fastq.gz
HPAP-171/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-171_GEX_NSS0032_5000904_S1_L002_R1_001.fastq.gz
HPAP-171/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-171_GEX_NSS0032_5000904_S1_L002_R2_001.fastq.gz
HPAP-172/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-172_GEX_NSS0032_5000906_S2_L001_I1_001.fastq.gz
HPAP-172/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-172_GEX_NSS0032_5000906_S2_L001_I2_001.fastq.gz
HPAP-172/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-172_GEX_NSS0032_5000906_S2_L001_R1_001.fastq.gz
HPAP-172/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-172_GEX_NSS0032_5000906_S2_L001_R2_001.fastq.gz
HPAP-172/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-172_GEX_NSS0032_5000906_S2_L002_I1_001.fastq.gz
HPAP-172/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-172_GEX_NSS0032_5000906_S2_L002_I2_001.fastq.gz
HPAP-172/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-172_GEX_NSS0032_5000906_S2_L002_R1_001.fastq.gz
HPAP-172/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-172_GEX_NSS0032_5000906_S2_L002_R2_001.fastq.gz
HPAP-173/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-173_GEX_NSS0032_5000908_S3_L001_I1_001.fastq.gz
HPAP-173/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-173_GEX_NSS0032_5000908_S3_L001_I2_001.fastq.gz
HPAP-173/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-173_GEX_NSS0032_5000908_S3_L001_R1_001.fastq.gz
HPAP-173/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-173_GEX_NSS0032_5000908_S3_L001_R2_001.fastq.gz
HPAP-173/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-173_GEX_NSS0032_5000908_S3_L002_I1_001.fastq.gz
HPAP-173/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-173_GEX_NSS0032_5000908_S3_L002_I2_001.fastq.gz
HPAP-173/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-173_GEX_NSS0032_5000908_S3_L002_R1_001.fastq.gz
HPAP-173/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-173_GEX_NSS0032_5000908_S3_L002_R2_001.fastq.gz
HPAP-174/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-174_GEX_NSS0032_5000910_S4_L001_I1_001.fastq.gz
HPAP-174/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-174_GEX_NSS0032_5000910_S4_L001_I2_001.fastq.gz
HPAP-174/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-174_GEX_NSS0032_5000910_S4_L001_R1_001.fastq.gz
HPAP-174/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-174_GEX_NSS0032_5000910_S4_L001_R2_001.fastq.gz
HPAP-174/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-174_GEX_NSS0032_5000910_S4_L002_I1_001.fastq.gz
HPAP-174/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-174_GEX_NSS0032_5000910_S4_L002_I2_001.fastq.gz
HPAP-174/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-174_GEX_NSS0032_5000910_S4_L002_R1_001.fastq.gz
HPAP-174/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-174_GEX_NSS0032_5000910_S4_L002_R2_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_GEX_NSS0032_5000912_S5_L001_I1_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_GEX_NSS0032_5000912_S5_L001_I2_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_GEX_NSS0032_5000912_S5_L001_R1_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_GEX_NSS0032_5000912_S5_L001_R2_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_GEX_NSS0032_5000912_S5_L002_I1_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_GEX_NSS0032_5000912_S5_L002_I2_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_GEX_NSS0032_5000912_S5_L002_R1_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_GEX_NSS0032_5000912_S5_L002_R2_001.fastq.gz
"

# Set IFS (Internal Field Separator)
IFS_BAK=$IFS
IFS=$'\n'

cd ./hpapdata

for f in $FILES; do
    echo "[$(date -Iseconds)] downloading $(basename $f) ..."
    encoded_f="$(echo $f | sed 's/ /%20/g')"
    curl --silent --create-dirs --output $f ${DATA_SERVER}/${encoded_f}
done

# Recover original IFS
IFS=$IFS_BAK
unset IFS_BAK

cd ..
echo; echo "[$(date -Iseconds)] experiment data downloaded"; echo

