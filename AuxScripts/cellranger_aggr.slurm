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

# This script uses a CSV file that is created in cellranger_count.slurm with sample IDs in first column and path to cellranger-count outputs in second column
# This CSV requires column names = "sample_id" and "molecule_h5" which are added at the beginning of this script

. ~/.profile
echo "cellranger-aggr"
echo ""

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME

source activate /global/smf/miniconda38_admin/envs/CellRanger_12-2021

glds=$1

outDir=/global/scratch/lmsande2/scRNAseq/GLDS-$glds/01-CellRanger/CellRanger_Output
# TODO: replace with permanent path to scRNAseq dir

fastqDir=/global/scratch/lmsande2/scRNAseq/GLDS-$glds/00-Raw_Data/Fastq
# TODO: replace with permanent path to scRNAseq dir

# Write out column headers and samples and their out dirs  to the aggr CSV files
# Samples and their out dirs
# grep: Grab only the R1 files since sample names are the same for R1 and R2 files
# awk: cut off everything after _S1 since everything before that is the sample name
for sample in $(ls $fastqDir | grep 'R1' | awk -F '_S1' '{print $1}'); do 
echo $sample','$outDir/$sample'/outs/molecule_info.h5' >> $outDir/GLDS-$glds\_aggr_CSV.csv
done
# Column headers
sed -i '1s/^/sample_id,molecule_h5\n/' $outDir/GLDS-$glds\_aggr_CSV.csv


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

