#!/bin/bash
#
#SBATCH --job-name="cellranger-mkref"
#SBATCH --output=cellranger-mkref.out
#SBATCH --partition=priority
#SBATCH --mem=50000
#SBATCH --mail-user=lauren.m.sanders@nasa.gov
#SBATCH --mail-type=END

source activate /global/smf/miniconda38_admin/envs/CellRanger_12-2021

cellranger mkref --genome=mmusculus_GRCm38.101_mkref \
		--fasta=/global/data/Data_Processing/Genomes/ensembl_release101/Mus_musculus/Mus_musculus.GRCm38.dna.primary_assembly.fa \
		--genes=/global/data/Data_Processing/Genomes/ensembl_release101/Mus_musculus/Mus_musculus.GRCm38.101.gtf

