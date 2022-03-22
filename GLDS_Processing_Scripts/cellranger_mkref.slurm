#!/bin/bash
#
#SBATCH --job-name="cellranger-mkref"
#SBATCH --output=cellranger-mkref.out
#SBATCH --partition=priority
#SBATCH --mem=50000
#SBATCH --mail-user=lauren.m.sanders@nasa.gov
#SBATCH --mail-type=END

. ~/.profile
echo "cellranger-mkref"
echo ""

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID


source activate /global/smf/miniconda38_admin/envs/CellRanger_12-2021

call="cellranger mkref --genome=mmusculus_GRCm38.101_mkref \
		--fasta=/global/data/Data_Processing/Genomes/ensembl_release101/Mus_musculus/Mus_musculus.GRCm38.dna.primary_assembly.fa \
		--genes=/global/data/Data_Processing/Genomes/ensembl_release101/Mus_musculus/Mus_musculus.GRCm38.101.gtf"


echo $call
eval $call
echo ""






echo ""
end=$(date +%s)	
echo "end time: $end"

echo ""
runtime=$(expr $end - $start)
echo "total runtime(s): $runtime"

sec_per_min=60
sec_per_hour=3600

runtime_m=$(echo "scale=2; $runtime / $sec_per_min;" | bc)
echo "total runtime(m): $runtime_m"

runtime_h=$(echo "scale=2; $runtime / $sec_per_hour;" | bc)
echo "total runtime(hr): $runtime_h"
