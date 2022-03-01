#!/bin/bash
#
#SBATCH --job-name="cellranger-count"
#SBATCH --output=cellranger-count.out
#SBATCH --partition=priority
#SBATCH --mem=50000
#SBATCH -c 80
#SBATCH --mail-user=lauren.m.sanders@nasa.gov
#SBATCH --mail-type=END


# Example script call: sbatch cellranger_count.slurm SampleID Path/To/Cellranger/Reference/Transcriptome /Path/To/Fastqs 8000 78 Path/To/Outdir
# This script writes the cellranger count outputs to a sample-named dir in the working directory, then moves it to the desired Outdir

. ~/.profile
echo "cellranger-count"
echo ""

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME

source activate /global/smf/miniconda38_admin/envs/CellRanger_12-2021

sample=$1
reference_transcriptome=$2
fastq_path=$3
expect_cells=$4
local_cores=$5
output_directory=$6

call="cellranger count --id=$sample \
--transcriptome=$reference_transcriptome \
--fastqs=$fastq_path \
--sample=$sample \
--expect-cells=$expect_cells \
--localcores=$local_cores 
"


echo $call
eval $call
echo “”

call1 = "mv $sample $output_directory"
eval $call1


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

