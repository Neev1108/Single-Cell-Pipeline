#!/bin/bash
#
#SBATCH --job-name="cellranger-aggr"
#SBATCH --output=cellranger-aggr.out
#SBATCH --partition=priority
#SBATCH --mem=50000
#SBATCH -c 80
#SBATCH --mail-user=lauren.m.sanders@nasa.gov
#SBATCH --mail-type=END


# Example script call: sbatch cellranger_aggr.slurm GLDS#
# This script runs cellranger-aggr to aggregate all the outputs and writes the filtered_feature_bc_matrix to 01-CellRanger/CellRanger_Output/GLDS#-A>

. ~/.profile
echo "cellranger-aggr"
echo ""

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME

source activate /global/smf/miniconda38_admin/envs/CellRanger_12-2021

glds=$1

outDir=/global/data/temp_scratch/lmsande2/scRNAseq/GLDS-$glds/01-CellRanger/CellRanger_Output
# TODO: replace with permanent path to scRNAseq dir

fastqDir=/global/data/temp_scratch/lmsande2/scRNAseq/GLDS-$glds/00-RawData/Fastq
# TODO: replace with permanent path to scRNAseq dir

# Run cellranger aggr
call="cellranger aggr --id=GLDS-$glds-AGGR \
--csv=$outDir/GLDS-$glds\_aggr_CSV.csv \
--nosecondary
"

echo $call
eval $call
echo ""

# Move the aggr output dir to the cellranger output dir
call2="mv GLDS-$glds-AGGR $outDir"
eval $call2

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

