## Scripts for downloading pancreas data

Scripts used to download pancreas data from [Pancreas DataBank](https://hpap.pmacs.upenn.edu/).

Files should be run in this order:
1. `create_genome_index.run`: generate STAR compatible genome index. can be skipped if you already have one
2. `pancdb_download_pt1.sh` and `pancdb_download_pt2.sh`: can be run to download `fastq.gz` files from PancDB
3. `generate_fastqc.sh` (optional) make fast-qc files to do quality control. Pancreas data is solid.
4. `run_star.sh` Used to run STARSolo on selected donors and create genome alignments
5. `generate_anndata.sh` Calls `generate_anndata_obj.py` to create a single anndata object from aligned output matrices of all of the selected donors
